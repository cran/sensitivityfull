senfmCI<-function(y,treated1,gamma=1,inner=0,trim=3,lambda=1/2,
                  alpha=0.05,twosided=TRUE,upper=TRUE){
  stopifnot(is.logical(treated1))
  stopifnot(gamma>=1)
  stopifnot((inner>=0)&(inner<=trim))
  stopifnot((lambda>0)&(lambda<1))
  stopifnot((alpha>0)&(alpha<1))
  stopifnot(is.logical(twosided))
  stopifnot(is.logical(upper))
  if(is.data.frame(y)) y<-as.matrix(y)
  stopifnot(is.matrix(y))
  stopifnot((dim(y)[1])==length(treated1))
  stopifnot(all(!is.na(as.vector(y[,1:2]))))

#Determine a finite interval, (-mx,mx) that will be searched for CI limits
  mx1<-max(y[,1],na.rm=TRUE)-min(y[,-1],na.rm=TRUE)
  mx2<-max(y[,-1],na.rm=TRUE)-min(y[,1],na.rm=TRUE)
  mx<-max(mx1,mx2)

#Internal function returning the tau that corresponds with a specified Normal quantile
  solve4tau<-function(dev,alternative="greater"){
    f<-function(taus){
      ntaus<-length(taus)
      o<-rep(NA,ntaus)
      for (i in 1:ntaus){
        d<-senfm(y,treated1,gamma=gamma,inner=inner,trim=trim,lambda=lambda,tau=taus[i],alternative=alternative)$deviate
        o[i]<-as.vector(d-dev)
      }
      o
    }
    rt<-uniroot(f,c(-mx,mx))
    unlist(rt$root)
  }

  if (twosided) nq<-(-qnorm(alpha/2))
  else nq<-(-qnorm(alpha))
  estL<-solve4tau(0,alternative="greater")
  if (gamma>1) estH<-solve4tau(0,alternative="less")
  else estH<-estL
  if (twosided|upper) ciL<-solve4tau(nq,alternative="greater")
  else ciL<-(-Inf)
  if (twosided|(!upper)) ciH<-solve4tau(nq,alternative="less")
  else ciH<-Inf

  CI<-"Lower"
  if(twosided) CI<-"Two-sided"
  if((!twosided)&upper) CI<-"Upper"
  desc<-c(1-alpha,gamma,CI)
  names(desc)<-c("Coverage","Gamma","Confidence Interval")

  list(PointEstimates=c(estL,estH),ConfidenceInterval=c(ciL,ciH),description=desc)
}
