library(sae)


#-------------------------------------------------#
#       Fitting transformed mixed model           #
#        with dual power transformation           #
#-------------------------------------------------#
### Input
# y: response vector
# X: matrix of sampled covariates (the first column should be 1's)
# rX: matrix of non-sampled covariates (the first column should be 1's)
# m: number of areas
# Ni: vector of the numbers of population in each area
# ni: vector of the numbers of samples in each area
# z: poverty line
# la.set: set of lambda
# mc: number of Monte Carlo samples
### Output
# Pred: estimated area-wise poverty measure
# Estimates: parameter estimates 
# RE: Estimates of random effects 
# MC: All Monte calro samples  
TMMP=function(y,X,rX,m,Ni,ni,z,la.set=seq(0,1,by=0.05),mc=1000){
  pov=function(u){ ifelse(u<z,1,0) }  # this can be customized
  p=dim(X)[2]; dom=rep(1:m,ni)
  mat1=cbind(1:m,matrix(1,m,p))
  mat2=cbind(1:m,ni)
  
  H=function(u,la){
    if(la==0){log(u)}
    else{ (u^la-u^(-la))/(2*la) }
  }
  invH=function(u,la){
    if(la==0){exp(u)}
    else{ (la*u+sqrt(la^2*u^2+1))^(1/la) }
  }
  
  plike=function(la){
    yy=H(y,la)
    res=eblupBHF(yy~X[,-1],dom=dom,meanxpop=mat1,popnsize=mat2,method="ML")
    Jac=log(y^(la-1)+y^(-la-1))-log(2)
    val=(res$fit$summary$logLik+sum(Jac))/sum(ni)
    return(val)
  }
  
  len=length(la.set)
  PL=c()
  for(k in 1:len){ PL[k]=plike(la.set[k]) }
  hla=la.set[which.max(PL)]
  
  yy=H(y,hla)
  res=eblupBHF(yy~X[,-1],dom=dom,meanxpop=mat1,popnsize=mat2,method="ML")
  hbeta=res$fit$summary$coefficient[,1]
  hsig=sqrt(res$fit$errorvar); htau=sqrt(res$fit$refvar)
  th=unlist(res$fit$random)
  hs=sqrt(hsig^2*htau^2/(hsig^2+ni*htau^2))
  
  if(hla==0){
    invH=function(u){ exp(u) }
  }else{
    invH=function(u){ (hla*u+sqrt(hla^2*u^2+1))^(1/hla) }
  }
  
  Est=c(hbeta,htau,hsig,hla)
  names(Est)=c(paste0("beta",1:p),"tau","sig","lambda")
  
  ri=Ni-ni
  rcu=c(0,cumsum(ri))
  scu=c(0,cumsum(ni))
  Ncu=c(0,cumsum(Ni))
  rmu=as.vector(rX%*%hbeta)
  MC=matrix(NA,mc,sum(Ni))
  Pred=c()
  
  for(i in 1:m){
    sterm=(scu[i]+1):scu[i+1]
    SY=t(matrix(y[sterm],ni[i],mc))
    
    posv=rnorm(mc,th[i],hs[i])
    rterm=(rcu[i]+1):rcu[i+1]
    RY=c()
    for(j in rterm){
      rn=rmu[j]+posv+rnorm(mc,0,hsig)
      RY=cbind(RY,invH(rn))
    }
    Rn=cbind(SY,RY)
    MC[,(Ncu[i]+1):Ncu[i+1]]=Rn
    Pred[i]=mean(pov(Rn))
  }
  
  Res=list(Pred,Est,th,MC)
  names(Res)=c("Predict","Estimates","RE","MC")
  return(Res)
}



#-------------------------------------------------#
#      Empirical Bayes confidence interval        #
#        (with dual power transformation)         #
#-------------------------------------------------#
### Input
# y: response vector
# X: matrix of sampled covariates (the first column should be 1's)
# rX: matrix of non-sampled covariates (the first column should be 1's)
# m: number of areas
# Ni: vector of the numbers of population in each area
# ni: vector of the numbers of samples in each area
# z: poverty line
# la.set: set of lambda
# mc: number of Monte Carlo samples
# alpha: significance level
# B: number of boostrap
### output
# area-wise confidence intervals
TMMP.CI=function(y,X,rX,m,Ni,ni,z,la.set=seq(0,1,by=0.05),mc=1000,alpha=0.05,B=100,print=T){
  pov=function(u){ ifelse(u<z,1,0) }  # this can be customized
  p=dim(X)[2]; dom=rep(1:m,ni)
  
  invH=function(u,la){
    if(la==0){exp(u)}
    else{ (la*u+sqrt(la^2*u^2+1))^(1/la) }
  }
  
  res=TMMP(y,X,rX,m,Ni,ni,z,la.set=la.set,mc=mc)
  est=res$Estimates; th=res$RE
  
  AM=function(x){
    am=c(); cu=c(0,cumsum(Ni))
    for(i in 1:m){ am[i]=mean(x[(cu[i]+1):cu[i+1]]) }
    return(am)
  }
  pos.mu=t(apply(pov(res$MC),1,AM))
  hbeta=est[1:p]; htau=est[p+1]; hsig=est[p+2]; hla=est[p+3]
  hs=sqrt(hsig^2*htau^2/(hsig^2+ni*htau^2))
  
  DGP.boot=function(para){
    beta=para[1:p]; tau=para[p+1]; sig=para[p+2]; la=para[p+3]
    mu=as.vector(X%*%beta); rmu=as.vector(rX%*%beta)
    cu=c(0,cumsum(ni))
    ri=Ni-ni; rcu=c(0,cumsum(ri))
    true=c(); y=c()
    v=rnorm(m,0,tau)
    for(i in 1:m){
      term=(cu[i]+1):cu[i+1]; rterm=(rcu[i]+1):rcu[i+1]
      rn1=mu[term]+v[i]+rnorm(ni[i],0,sig)
      rn2=rmu[rterm]+v[i]+rnorm(ri[i],0,sig)
      y[term]=invH(rn1,la); Y=invH(c(rn1,rn2),la)
      true[i]=mean(pov(Y))
    }
    return(list(y,true))
  }
  
  boot.sample=list(); boot.true=list()
  for(b in 1:B){
    boot=DGP.boot(est)
    boot.y=boot[[1]]; boot.true[[b]]=boot[[2]]
    boot.res=TMMP(boot.y,X,rX,m,Ni,ni,z,la.set=seq(hla-0.1,hla+0.1,by=0.05),mc=mc)
    boot.sample[[b]]=t(apply(pov(boot.res$MC),1,AM))
    if(print) { print(paste0("Bootstrap replication (b=",b,")")) }
  }
  
  len=11
  a.set=seq(0,alpha,length=len)
  cv.prob=c()
  for(k in 1:len){
    a=a.set[k]
    qq=function(x){ quantile(x,prob=c(a/2,1-a/2)) }
    
    id=matrix(NA,B,m)
    for(b in 1:B){
      CI=apply(boot.sample[[b]],2,qq)
      id[b,]=ifelse(CI[1,]<boot.true[[b]],1,0)*ifelse(CI[2,]>boot.true[[b]],1,0)
    }
    cv.prob[k]=mean(id)
  }
  
  opt.a=a.set[which.min(abs(cv.prob-(1-alpha)))]
  calpha=opt.a
  cquant=function(x){ quantile(x,prob=c(calpha/2,1-calpha/2)) }
  CI=apply(pos.mu,2,cquant)
  dimnames(CI)[[1]]=c(alpha/2,1-alpha/2)
  dimnames(CI)[[2]]=paste0("area-",1:m)
  return(t(CI))
}



