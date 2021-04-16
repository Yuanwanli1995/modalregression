# 众数回归模型设计
# 生成众数回归模型类
# 模型训练
datagenerater<- function(n,e,f){
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
  }else if(e=="mixgausssymmetry"){# 重尾
    epsilon1<-rnorm(0.9*n,mean=0,sd=1)
    epsilon2<-rnorm(0.05*n,mean=4,sd=1)
    epsilon3<-rnorm(0.05*n,mean=-4,sd=1)
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
  }else if(f=="linear"){
    func<-function(x){return(2*x+5)}#seed(1)
  }
  y<-func(x)+error
  out <- list(x=data.frame(x),y=data.frame(y),real=data.frame(func(x)+mode))
  return(out)
}
l2<-function(c){return(sqrt(t(c)%*%c))}# 计算向量范数
gausskernel<-function(t1,t2,h){
  return(exp((t1-t2)^2/(-2*h^2))/(sqrt(2*pi)*h))}# 高斯核函数
kernel_dm<-function(x,x_data){# 生成lmr设计矩阵
  n<-length(x_data);
  d<-length(x)
  dm<-matrix(data = rep(0,n*d), nrow = n);
  for(i in seq(n)){
    for(j in seq(d)){
      dm[i,j]<-gausskernel(x_data[i],x[j],0.1)
    }
  }
  return(dm)
}
 
ite<-function(dm,y,beta,sp_1,sp_2){
  #beta ->new beta  # liner,bmr,kmr迭代函数
  dimension<-ncol(dm)
  response<-gausskernel(y,dm%*%beta,sp_1)
  W <- diag(as.vector(response/sum(response)))
  temp<-t(dm)%*%W%*%dm+diag(rep(sp_2,dimension))
  x.inv <- try(solve(temp),silent=T)
  if ('try-error' %in% class(x.inv)) return(beta)
  else return(solve(temp)%*%t(dm)%*%W%*%y)
}
ite_lmr<-function(method,x,y,dm,beta,sp_1,sp_2,x0){
  dimension<-ncol(dm)
  response<-gausskernel(y,dm%*%beta,sp_1)*gausskernel(x,x0,sp_2)
  W <- diag(as.vector(response/sum(response)))
  temp<-t(dm)%*%W%*%dm
  x.inv <- try(solve(temp),silent=T)
  if ('try-error' %in% class(x.inv)) return(beta)
  else return(solve(temp)%*%t(dm)%*%W%*%y)
}
setClass(
  "MR",
  slots=list(method="character",
             x_train="data.frame",
             y_train="data.frame",
             hyperparameters="numeric",
             n="numeric",
             d="numeric",
             coefficients="numeric",
             fit="numeric",
             iterations="numeric",
             bandwidth_criterion="character",
             interval1="numeric",
             interval2="numeric",
             mae="numeric",
             computing_time="numeric"),
  prototype =list(method="linear",
                  hyperparameters=NaN,
                  coefficients=NaN,
                  fit=NaN,
                  iterations=0,
                  bandwidth_criterion="mse",
                  interval1=seq(0.05,1,0.05),
                  interval2=c(0.1,1,10,100,1000),
                  mae=NaN,
                  computing_time=NaN)
)
setValidity("MR", function(object){
  if (!object@method %in% c("linear","lmr","bmr","kmr"))
    stop("Sorry,  The method you choose is not currently supported,
    linear、lmr、bmr or kmr are available")
  if(!any(is.nan(object@hyperparameters))){
    if(any(object@hyperparameters<0))
      stop("hyperparameters should be positive")}
  if (!object@bandwidth_criterion %in% c("mse","predict"))
    stop("Sorry,  The bandwidth selection criterion you choose is not currently supported,
    mse、predict are available")
})
# 模型实现
# 最多输入四个参数， 最少两个
# 输出一个mr模型类
modalreg<-function(x_train, y_train, method="bmr", 
                   bandwidth_criterion="mse", hyperparameters=NaN,
                   interval1=seq(0.05,1,0.05),interval2=c(0.1,1,10,100,1000)){
  ptm<- proc.time()[3]# 开始时间
  n = nrow(x_train);
  d = ncol(x_train);
  # 实例化一个MR对象
  mrobject<-new("MR",
                method=method,
                x_train=x_train,
                y_train=y_train,
                hyperparameters = hyperparameters,
                n=n,
                d=d,
                bandwidth_criterion=bandwidth_criterion,
                interval1=interval1,
                interval2=interval2)
  # 如果没有提供超参数，自动选择； 否则，使用. 但为负数检查不通过。
  if(any(is.nan(mrobject@hyperparameters))){
    # bandwidth函数面向mr， 返回最优超参数向量
    mrobject@hyperparameters<-bandwidth(mrobject)
  }
  # mrfit函数面向mr，返回回归系数,和训练集预测结果,迭代次数
  out<-mrfit(mrobject)
  mrobject@coefficients<-out$coefficients
  mrobject@fit<-out$fit
  mrobject@mae<-mean(abs(as.matrix(y_train)-out$fit))
  mrobject@iterations<-out$iterations
  mrobject@computing_time=proc.time()[3]-ptm
  return(mrobject)
}
# 训练模型
# 设定模型拟合泛型函数
setGeneric("mrfit",
           function(object){
             standardGeneric("mrfit")
           }
)
### mrfit: data, papa-> beta, ceo, fit 
setMethod("mrfit",
          signature(object = "MR"),
          function(object){
            # 解包
            x<-as.matrix(object@x_train);
            y<-as.matrix(object@y_train);
            method<-object@method
            sp_1<-object@hyperparameters[1]#sp1为kde核，
            sp_2<-object@hyperparameters[2]#sp2为光滑参数
            n<-object@n
            d<-object@d
            # 回归
            if(method=="linear"){
              # design matrix 
              dm<-cbind(rep(1,n),x)
              beta = matrix(rep(0,ncol(dm)))
              iterations <- 0 
              repeat{
                beta_next<-ite(dm,y,beta,sp_1,0)
                e=l2(beta_next-beta)
                if(e>0.000001){
                  beta<-beta_next
                  iterations = iterations+1
                }else{break}
              }
              coefficients=as.vector(beta);fit=as.vector(dm%*%beta)
            }else if(method=="bmr"){
              if(d!=1){stop("非参数方法仅支持一维")}
              library(splines)
              aa=seq(5,n-1,length.out = sp_2)
              knots1=c(0,0,0,0,x[aa],1,1,1,1) 
              dm<-splineDesign(knots=knots1, x, outer.ok = T)
              beta = matrix(rep(0,ncol(dm)))
              iterations <- 0 
              repeat{
                beta_next<-ite(dm,y,beta,sp_1,0)
                e=l2(beta_next-beta)
                if(e>0.000001){
                  beta<-beta_next
                  iterations = iterations+1
                }else{break}
              }
              coefficients=as.vector(beta);fit=as.vector(dm%*%beta)
            }else if(method =="kmr"){
              if(d!=1){stop("非参数方法仅支持一维")}
              dm<-kernel_dm(x,x)
              beta = matrix(rep(0,ncol(dm)))
              iterations <- 0 
              repeat{
                beta_next<-ite(dm,y,beta,sp_1,sp_2)
                e=l2(beta_next-beta)
                if(e>0.000001){
                  beta<-beta_next
                  iterations = iterations+1
                }else{break}
              }
              coefficients=as.vector(beta);fit=as.vector(dm%*%beta)
            }else if(method=="lmr"){
              if(d!=1){stop("非参数方法仅支持一维")}
              fit<-rep(0,n)
              iterations <- 0 
              for (i in seq(1,n)){
                x0<-x[i]
                dm<-cbind(rep(1,n),x-x0,(x-x0)^2)
                beta = matrix(rep(0,ncol(dm)))
                repeat{
                  beta_next<-ite_lmr(method,x,y,dm,beta,sp_1,sp_2,x0)
                  e=l2(beta_next-beta)
                  if(e>0.000001){
                    beta<-beta_next
                    iterations = iterations+1
                  }else{break}
                }
                fit[i]<-beta[1]
              }
              coefficients=NaN;
            }
            out <- list(coefficients=coefficients,fit=fit,iterations=iterations)
            return(out)
          }
)
setGeneric("bandwidth",
           function(object){
             standardGeneric("bandwidth")
           }
)
setMethod("bandwidth",
          signature(object="MR"),
          function(object){
            # 解包
            interval1<-object@interval1;
            interval2<-object@interval2;
            x<-as.matrix(object@x_train);
            y<-as.matrix(object@y_train);
            method<-object@method
            bandwidth_criterion<-object@bandwidth_criterion
            #选参数
            c<-(max(y)-min(y))*0.05
            S<<-matrix(rep(0,length(interval1)*length(interval2)),nrow = length(interval1))
            for(i in 1:length(interval1)){
              for(j in 1:length(interval2)){
                hyperparameters<-c(interval1[i],interval2[j])
                object@hyperparameters<-hyperparameters
                out<-mrfit(object)
                fit<-out$fit
                if(bandwidth_criterion=="predict"){
                  S[i,j]<<-sum(fit-c<y&y<fit+c)
                }else if(bandwidth_criterion=="mse"){
                  S[i,j]<<-sum((fit-y)^2)*(-1)
                }
              }
            }
            b<-as.matrix(which(S==max(S),arr.ind=T))
            sp_1<-interval1[b[1,1]]
            sp_2<-interval2[b[1,2]]
            return(c(sp_1,sp_2))
          })
modalpredict<-function(model, x_predict){
  if( class(x_train)!="data.frame"){
    stop("x_predict 类型需为数据框")
  }
  if( class(model)!="MR"){
    stop("model 类型需为MR类")
  }
  if( model@d!=ncol(x_predict)){
    stop("预测数据需要与训练数据维度一致")
  }
  x<-as.matrix(x_predict)
  x_train<-as.matrix(model@x_train)
  y_train<-as.matrix(model@y_train)
  n <-nrow(x_predict)
  d <-ncol(x_predict)
  method <-model@method
  beta <- model@coefficients
  sp_1 <- model@hyperparameters[1]
  sp_2 <- model@hyperparameters[2]
  if(method=="linear"){
    # design matrix
    dm<-cbind(rep(1,n),x)
    fit<-dm%*%beta 
  }else if(method=="bmr"){
    aa=seq(5,n-1,length.out = sp_2)
    knots1=c(0,0,0,0,x_train[aa],1,1,1,1)
    dm<-splineDesign(knots=knots1, x, outer.ok = T)
    fit<-dm%*%beta
  }else if(method =="kmr"){
    dm<-kernel_dm(x_train,x)
    fit<-dm%*%beta
  }else if(method=="lmr"){
    fit<-rep(0,n)
    for (i in seq(1,n)){
      x0<-x[i]
      dm<-cbind(rep(1,n),x_train-x0,(x_train-x0)^2)
      beta = matrix(rep(0,ncol(dm)))
      repeat{
        beta_next<-ite_lmr(method,x_train,y_train,dm,beta,sp_1,sp_2,x0)
        e=l2(beta_next-beta)
        if(e>0.000001){
          beta<-beta_next
        }else{break}
      }
      fit[i]<-beta[1]
    }
    coefficients=NaN;
  }
  return(fit)
}
setGeneric("modal.summary",
           function(object){
             standardGeneric("modal.summary")
           }
)
setMethod("modal.summary",
          signature(object = "MR"),
          function(object){
            if(object@iterations==0){
              cat("参数选择不合适，请尝试调整参数或搜索区间")
            }
            cat("n:", object@n)
            cat("          d:", object@d,"\n")
            cat("method:", object@method,"\n")
            cat("hyperparameters:", object@hyperparameters,"\n")
            if(object@method!="kmr"){
              cat("coefficients:", object@coefficients,"\n")
            }else{cat("coefficients:", object@coefficients[1:5],"  ...\n")}
            cat("iterations:", object@iterations,"\n")
            cat("mean abslute error:",object@mae,"\n")
            cat("computing time:",object@computing_time,"s","\n")
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
            if(object@iterations>0){
              y<-object@y_train;
              y<-y[[names(y)]]
              plot(seq(object@n),y,xlab="x",ylab="y_train")
              lines(seq(object@n),object@fit)
            }else{print("参数选择不合适，请尝试调整参数或搜索区间")}
          }
)
########example
# data<- datagenerater(200,"gamma","exp")
# x_train<-data$x
# y_train<-data$y
# bmr1<-modalreg(x_train = x_train, y_train = y_train, method="bmr", hyperparameters=c(2,1))
# bmr2<-modalreg(x_train = x_train, y_train = y_train, method="bmr",bandwidth_criterion="predict")
# bmr3<-modalreg(x_train = x_train, y_train = y_train, method="bmr",bandwidth_criterion="predict",
#                interval1=seq(1,10,1), interval2 = seq(0,1,0.25))
# kmr1<-modalreg(x_train = x_train, y_train = y_train, method="kmr",bandwidth_criterion="mse",
#                interval1=seq(1,10,1), interval2 = seq(0,1,0.5))
# lmr1<-modalreg(x_train = x_train, y_train = y_train, method="lmr",bandwidth_criterion="mse")
# modal.summary(bmr1)
# modal.summary(bmr2)
# modal.summary(bmr3)
# modal.summary(kmr1)
# modal.summary(lmr1)
# modal.plot(bmr1)
# pred<-modalpredict(model, x_train)
# plot(as.matrix(x_train), as.matrix(y_train),xlab ="x", ylab = "y")
# lines(as.matrix(x_train),pred,col="red")
