# 众数回归模型设计
# 生成众数回归模型类
library(splines)
generatedata<- function(n,e,f){
  n<-n
  x<-sort(runif(n,0,1))
  if(e=="gauss"){
    error <- rnorm(n, mean= 0, sd = 1)
    mode<-0
  }else if(e=="mixgauss"){
    epsilon1<-rnorm(0.1*n,mean=4,sd=1)
    epsilon2<-rnorm(0.9*n,mean=0,sd=1)
    epsilon<-c(epsilon1,epsilon2)
    error<-sample(epsilon,n,replace=TRUE)#biased error
    mode<-0
  }else if(e=="mixgausstwoside"){# 重尾
    epsilon1<-rnorm(0.9*n,mean=0,sd=1)
    epsilon2<-rnorm(0.05*n,mean=5,sd=1)
    epsilon3<-rnorm(0.05*n,mean=-5,sd=1)
    epsilon<-c(epsilon1,epsilon2,epsilon3)
    error<-sample(epsilon,n,replace=TRUE)#biased error
    mode<-0
  }else if(e=="cauchy"){ # 离群点
    error <- rcauchy(n, location = 0, scale = .1)
    mode<-0
  }else if(e == "gamma"){# 偏态
    error <- rgamma(n, shape=1,rate = 1)-0.3 # mode = 0 ,mean =0.7 ,sd =1
    mode<-0
  }else if(e == "beta"){# 可以是有偏
    error <- rbeta(n,5,3)
    mode<-0.7
  }
  if(f=="quodratic" ){
    func<-function(x){return(40*x^2-40*x+10)}##seed(6)
  }else if(f=="exp"){
    func<-function(x){return(x+exp(8*(x-0.5)^2))}
  }else if(f=="log"){
    func<-function(x){return(5*log(x,3)+10)}###logarithm
  }else if(f=="sin"){
    func<-function(x){return(5*sin(20*x))}#seed(43)
  }else if(f=="mexicohat"){
    func<-function(x){return(-1+1.5*x+0.2*dnorm(x-0.6,0,0.04))}#seed(1)
  }
  y<-func(x)+error
  out <- list(x=x,y=y,real=func(x)+mode)
  return(out)
}
# 模型训练
l2<-function(c){return(sqrt(t(c)%*%c))}#functions for EM
phih<-function(i,X,y,beta_bmr,h){
  return(exp((y[i]-X[i,]%*%beta_bmr)^2/(-2*h^2))/(sqrt(2*pi)*h))}
ite_bmr<-function(x,X,y,beta_bmr,h){
  n<-length(x)
  W<-as.vector(phih(1:n,X=X,y=y,beta_bmr=beta_bmr,h=h))
  W<-diag(W)
  x.inv <- try(solve(t(X)%*%W%*%X),silent=T)
  if ('try-error' %in% class(x.inv)) return(beta_bmr)
  else return(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y)}
lmr<-function(x,y,beta_lmr,x0,z,h1,h2){
  pre<-beta_lmr[1]*z+beta_lmr[2]*(x-x0)+beta_lmr[3]*(x-x0)^2
  return(exp(-((x-x0)^2)/(2*h1^2)-((y-pre)^2)/(2*h2^2)))}
ite_lmr<-function(x,X,y,beta_lmr,x0,z,h1,h2){
  W<-diag(lmr(x,y,beta_lmr,x0,z,h1,h2))
  x.inv <- try(solve(t(X)%*%W%*%X),silent=T)
  if ('try-error' %in% class(x.inv)) return(beta_lmr)
  else return(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y)}
regression<-function(x,y,method,sp_1, sp_2){
  if(method=="bmr"){
    h<-sp_1
    a<-sp_2
    aa=seq(0,n-1,a)
    knots1=c(0,0,0,0,x[aa],1,1,1,1)
    X<-splineDesign(knots=knots1,x, outer.ok = F)
    beta_bmr<-rep(0,ncol(X))
    repeat{
      temp<-ite_bmr(x,X,y,beta_bmr,h)
      e=l2(temp-beta_bmr)
      beta_bmr<-temp
      if(e<0.001)break}
    fit_bmr<-X%*%beta_bmr
    fit<-fit_bmr
  }else if(method=="lmr"){
    fit_lmr<-vector()
    z<-rep(1,n)
    h1<-sp_1
    h2<-sp_2
    beta_lmr<-rep(1,3)
    for(k in 1:n){
      x0<-x[k]
      X<-cbind(z,x-x0,(x-x0)^2)
      repeat{
        temp<-ite_lmr(x,X,y,beta_lmr,x0,z,h1,h2)
        e=l2(temp-beta_lmr)
        beta_lmr<-temp
        if(e<0.01)break}
      fit_lmr[k]<-beta_lmr[1]}
    fit<- fit_lmr
  }else if(method=="bqr"){
    fit<-fitted(rq(y~bs(x,df=3),tau=.5))
  }else if(method=="blse"){
    fit<-fitted(lm(y~bs(x,df=3)))
  }
  return(fit)
}
bandwidthselecte<-function(x,y,model,criterion ){
  c<-(max(y)-min(y))*0.05
  if(model=="bmr"){
    A<-seq(20,80,by=10)
    H<-seq(0.1,2,by=0.2)
    S<<-matrix(rep(0,length(A)*length(H)),nrow = length(H))
    for(i in 1:length(H)){
      h<-H[i]
      for(j in 1:length(A)){
        a<-A[j]
        fit_bmr<-regression(x,y,method="bmr",h,a)
        if(criterion=="p"){
          S[i,j]<-sum(fit_bmr-c<y&y<fit_bmr+c)
        }else if(criterion=="mse"){
          S[i,j]<-sum((fit_bmr-y)^2)*(-1)
        }
      }}
    #print(S)
    b<-as.matrix(which(S==max(S),arr.ind=T))
    sp_1<-H[b[1,1]]
    sp_2<-A[b[1,2]]
  }else if(model=="lmr"){
    H1<-seq(0.05,0.2,by=0.02)#h2 search and h1 prety small[0.01,0.1]
    H2<-seq(0.01,1,by=0.1)#h2 search
    S1<-matrix(rep(0,length(H1)*length(H2)),nrow = length(H1))
    for(i in 1:length(H1)){
      h1<-H1[i]
      for(j in 1:length(H2)){h2<-H2[j]
      fit_lmr<-regression(x,y,method="lmr",h1,h2)
      if(criterion=="p"){
        S1[i,j]<-sum(fit_lmr-c<y&y<fit_lmr+c)
      }else if(criterion=="mse"){
        S1[i,j]<-sum((fit_lmr-y)^2)*(-1)
      }
      }}
    # print(S1)
    b1<-as.matrix(which(S1==max(S1),arr.ind=T))
    sp_1<-H1[b1[nrow(b1),1]]
    sp_2<-H2[b1[nrow(b1),2]]
  }else if(model=="blse"){
    sp_1<--1;sp_2<--1
  }else if(model=="bqr"){
    sp_1<--2;sp_2<--2
  }
  out <- list(sp_1=sp_1,sp_2=sp_2)
  return(out)
}
setClass(
  "MR",
  slots=list(method="character",
             x_train="data.frame",
             y_train="data.frame",
             hyperparameters="numeric",
             n="numeric",
             d="numeric",
             Coefficients="numeric",
             fit="numeric"),
  prototype =list(method="linear",
                  hyperparameters=c(1,1),
                  Coefficients=c(0),
                  fit=c(0))
)
setValidity("MR", function(object){
  if (!object@method %in% c("linear","lpm","bmr","kmr"))
    stop("Sorry,  The method you choose is not currently supported,
    linear、lpm、bmr or kmr available")
})
# 模型实现
# 最多输入四个参数， 最少两个
# 输出一个mr模型类
modalreg<-function(x_train, y_train, method="linear", criterion="lse",hyperparameters=c(1,1)){
  n = nrow(x_train);
  d = ncol(x_train);
  # 实例化一个MR对象
  mrobject<-new("MR",
                 method=method,
                 x_train=x_train,
                 y_train=y_train,
                 hyperparameters = hyperparameters,
                 n=n,d=d)
  # 如果提供超参数，自动选择
  if(any(hyperparameters!= c(0.5,50))){
    # bandwidth函数面向mr， 返回最优超参数向量
    mrobject@hyperparameters<-bandwidth(mrobject)
  }
  # mrfit函数面向mr，返回回归系数,和训练集预测结果
  out<-mrfit(mrobject)
  mrobject@Coefficients<-out$para
  print("ok2")
  mrobject@fit<-out$fit
  return(mrobject)
}

# 训练模型
# 设定模型拟合泛型函数
setGeneric("mrfit",
           function(object){
             standardGeneric("mrfit")
           }
)
setMethod("mrfit",
          signature(object = "MR"),
          function(object){   #//function(x,y,method,sp_1, sp_2){//
            # 解包
            x<-object@x_train;
            x<-x[[names(x)]]
            y<-object@y_train;
            y<-y[[names(y)]]
            method<-object@method
            sp_1<-object@hyperparameters[1]
            sp_2<-object@hyperparameters[2]
            n<-object@n
            # 回归
            if(method=="bmr"){
              h<-sp_1
              a<-sp_2
              aa=seq(0,n-1,a)
              knots1=c(0,0,0,0,x[aa],1,1,1,1)
              X<-splineDesign(knots=knots1,x, outer.ok = F)
              beta_bmr<-rep(0,ncol(X))
              repeat{
                temp<-ite_bmr(x,X,y,beta_bmr,h)
                e=l2(temp-beta_bmr)
                beta_bmr<-temp
                if(e<0.001)break}
              para<-beta_bmr
              fit_bmr<-X%*%beta_bmr
              fit<-fit_bmr
            }else if(method=="lmr"){
              fit_lmr<-vector()
              z<-rep(1,n)
              h1<-sp_1
              h2<-sp_2
              beta_lmr<-rep(1,3)
              for(k in 1:n){
                x0<-x[k]
                X<-cbind(z,x-x0,(x-x0)^2)
                repeat{
                  temp<-ite_lmr(x,X,y,beta_lmr,x0,z,h1,h2)
                  e=l2(temp-beta_lmr)
                  beta_lmr<-temp
                  if(e<0.01)break}
                fit_lmr[k]<-beta_lmr[1]}
              para<-beta_lmr
              fit<- fit_lmr
            }
            print("ok")
            print(class(as.vector(fit)))
            print(class(as.vector(para)))
            out <- list(fit=as.vector(fit), para=as.vector(para))
            return(out)
          }
)
setGeneric("modal.summary",
           function(object){
             standardGeneric("modal.summary")
           }
)
setMethod("modal.summary",
          signature(object = "MR"),
          function(object){
            cat("method:", object@method,"\n")
            cat("hyperparameters:", object@hyperparameters,"\n")
            cat("Coefficients:", object@Coefficients,"\n")
          }
)
setGeneric("modal.plot",
           function(object){
             standardGeneric("modal.plot")
           }
)
setMethod("modal.plot",
          signature(object = "MR"),
          function(object){
            y<-object@y_train;
            y<-y[[names(y)]]
            print(y)
            plot(seq(object@n),y)
            lines(seq(object@n),object@fit)
          }
)
########example
#data<-generatedata(200,"gauss","sin")
#model<-modalreg(x_train = data.frame(data$x),y_train = data.frame(data$y),
#         method = "bmr",hyperparameters= c(0.5,50))
#modal.summary(model)
#modal.plot(model)
