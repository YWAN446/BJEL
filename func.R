#functions;
m_linear <- function(X, Epsilon){
  Y <- 1 + X + Epsilon
  return(Y)
}

m_nonlinear <- function(X, Epsilon){
  Y <- 1 + X^2 + Epsilon
  return(Y)
}

logit_linear <- function(X){
  pi <- 1/(1+exp(-1-X))
  return(pi)
}

logit_nonlinear <- function(X){
  pi <- 1/(1+exp(-1-X^2))
  return(pi)
}

Lag=function(u,ds,mu){
  dif=1
  tol=1e-8
  if(min(u-mu)>=0 | max(u-mu)<=0){
    dif=0
    M=0
  }
  L=-1/max(u-mu)
  R=-1/min(u-mu)
  while(dif>tol){
    M=(L+R)/2
    glam=sum((ds*(u-mu))/(1+M*(u-mu)))
    if(glam>0) L=M
    if(glam<0) R=M
    dif=abs(glam)
    #cat(paste0(M," ",dif,"\n"))
  }
  return(M)
}

#function to solve the Langrange multiplier;
Lfs2<-function(u,ds){
  dss<-ds/sum(ds)
  lambda0<-0
  k<-0
  gamma0<-1
  epsilon<-10^{-10}
  nu<-length(u)
  
  repeat{
    D1<-sum(dss*u/(1+lambda0*u))
    D2<-D1/(-sum(u^{2}*dss/(1+lambda0*u)^{2}))
    if(abs(D2)<epsilon){M<-lambda0
    break
    }
    repeat{ 
      deta0<-gamma0*D2
      if(sum(as.numeric(1+(lambda0-deta0)*u<=0))>0){ gamma0<-gamma0/2
      }
      else{break}
    }
    lambda0<-lambda0-deta0
    k<-k+1
    gamma0<-(k+1)^{-1/2}
    if(k>=100){ M<-0
    break
    }
  }
  return(M)
}

#####4.Lambda function for vector:

###u is a r*n matrix 
Lfv<-function(u){
  
  lambda0<-rep(0,dim(u)[1])
  k<-0
  gamma0<-1
  epsilon<-10^{-8}
  nu<-dim(u)[2]
  
  repeat{
    a<-t(1+t(lambda0)%*%u)
    a<-as.numeric(a)
    B<-NULL
    iter<-0
    repeat{
      iter<-iter+1
      B<-rbind(B,a)
      if(iter==dim(u)[1]){break}
    }
    unew<-u/B
    D1<-apply(unew,1,mean)
    D2<-solve(-nu^{-1}*unew%*%t(unew))%*%D1
    if(sqrt(sum(D2^{2}))<epsilon){M<-lambda0
    break
    }
    repeat{ 
      deta0<-gamma0*D2
      if(sum(as.numeric(t(1+t(lambda0-deta0)%*%u)<=0))>0){ gamma0<-gamma0/2
      }
      else{break}
    }
    lambda0<-lambda0-deta0
    k<-k+1
    ###    gamma0<-(k+1)^{-1/2}
    gamma0<-1
    if(k>=200){ M<-rep(0,dim(u)[1])
    break
    }
  }
  M<-as.numeric(M)
  return(M)
}

#generate the data;
##################################################
#input parameters:################################
##B:Monte Carlo sample size#######################
##N:sample size########################
#output results:##################################
##3 dimensional matrix contains y,x,delta#########
##################################################
gedata<-function(B,N,n,seed){
  
  ###super population model:
  set.seed(seed[1])
  x1<-rchisq(B*N,1)-1
  x2<-rchisq(B*N,3)-3
  x3<-rnorm(B*N,0,1)
  
  set.seed(seed[2])
  #epsilon<-rnorm(B*N,0,1)
  epsilon<-rchisq(B*N,1)-1
  
  y <- m_linear(x1+x2, epsilon)
  
  p <- logit_linear(x1+x2)
  
  set.seed(seed[3])
  delta<-rbinom(B*N,1,p) #delta==1 is respondent and delta==0 is not respondent;
  
  sI <- c()
  for (i in 1:B){
    temp <- rep(0,N)
    temp[sample(1:N,n)] <- 1
    sI <- c(sI,temp)
  }
  
  My<-matrix(y,B,N,byrow=T)
  Mx1<-matrix(x1,B,N,byrow=T)
  Mx2<-matrix(x2,B,N,byrow=T)
  Mx3<-matrix(x3,B,N,byrow=T)
  Mdelta<-matrix(delta,B,N,byrow=T)
  MsI<-matrix(sI,B,N,byrow=T)
  
  M<-cbind(My,Mx1,Mx2,Mx3,Mdelta,MsI)
  M<-t(M)
  aM<-as.numeric(M)
  Res<-array(aM,c(N,6,B))
  theta0<-c(mean(y),median(y))
  res<-list(Res,theta0)
  return(res)
}

#############functions for point estimation:

###4. Semiparametric fractional imputation estimator by using sample empirical likelihood:
fSFIS<-function(dat){
  num.col <- length(dat[1,])
  y<-dat[,1]
  x<-dat[,2:(num.col-2)] #dimension change
  delta<-dat[,num.col-1]
  sI<-dat[,num.col]
  n<-sum(sI)
  
  sx<-x[sI==1,] #dimension change
  sy<-y[sI==1]
  sdelta<-delta[sI==1]
  sN<-length(sy)
  
  rx<-sx[sdelta==1,] #dimension change
  ry<-sy[sdelta==1]
  nr<-length(which(sdelta==1))
  rw<-rep(1,nr)
  #rw<-rep(1,length(rx[,1])) #dimension change
  #nr<-length(rx[,1]) #dimension change
  rN<-length(ry)
  
  mx<-sx[sdelta==0,] #dimension change
  nm<-length(which(sdelta==0))
  mw<-rep(1,nm)
  #mw<-rep(1,length(mx[,1])) #dimension change
  #nm<-length(mx[,1]) #dimension change
  
  ###compute weighted least square estimator of beta:
  glm.dat <- data.frame(cbind(ry,rx))
  names(glm.dat) <- c("ry","X1","X2")
  model <- glm(formula = ry ~ ., data=glm.dat)
  em<-predict(model,data.frame(sx),type="response") #dimension change
  mem<-predict(model,data.frame(matrix(mx,ncol = 2)),type="response") #dimension change
  rem<-predict(model,data.frame(matrix(rx,ncol = 2)),type="response") #dimension change
  
  ###compute EL weights:
  
  re<-ry-rem
  u1<-rw*re
  elambda<-Lfv(t(u1))
  #elambda<-Lag(u1,rep(1,length(u1)),0)
  reg<-1/nr*1/(1+elambda*u1)
  reg2<-reg*rw/sum(reg*rw)
  
  ###compute the point estimate:
  
  ##matrix for eyij:
  a1<-rep(mem,nr)
  A1<-matrix(a1,nm,nr)
  b1<-rep(re,nm)
  B1<-matrix(b1,nm,nr,byrow=T)
  EYIJ<-A1+B1
  IEYIJ<-matrix(as.numeric(EYIJ<1),nm,nr)
  wg<-mw%o%reg2
  
  #mean
  etheta1<-(sum(rw*ry)+sum(EYIJ*wg))/sN
  
  #median
  kcdf <- function(x) {
    temp.x <- 0.75*x-0.25*x^3+0.5
    temp.x[x>=1] <- 1
    temp.x[x<=-1] <- 0
    return(temp.x)
  }
  
  #h <- h.amise(ry, deriv.order=0, kernel="epanechnikov")[[6]]
  h <- 1.06*sd(ry)/n^(1/5)
  fq<-function(tt){
    IEYIJ3<-matrix(kcdf((EYIJ-tt)/h),nm,nr)
    res<-((sum(rw*kcdf((ry-tt)/h))+sum(IEYIJ3*wg))/sN-0.5)^{2}
    return(res)
  }
  etheta2<-optimize(fq,quantile(y,c(0.05,0.95)))$minimum
  # if (abs(optimize(fq,c(min(y),max(y)))$minimum-max(y))<0.1){
  #   etheta2<-optimize(fq,c(min(y),theta0[2]+10))$minimum
  # } else {
  #   etheta2<-optimize(fq,c(min(y),max(y)))$minimum
  # }
  
  return(c(etheta1,etheta2))
}

###5. Propensity score weighted estimator:
fPS<-function(dat){
  num.col <- length(dat[1,])
  y<-dat[,1]
  x<-dat[,2:(num.col-2)] #dimension change
  delta<-dat[,num.col-1]
  sI<-dat[,num.col]
  n<-sum(sI)
  
  sx<-x[sI==1,] #dimension change
  sy<-y[sI==1]
  sdelta<-delta[sI==1]
  sN<-length(sy)
  
  rx<-sx[sdelta==1,] #dimension change
  ry<-sy[sdelta==1]
  nr<-length(which(sdelta==1))
  rw<-rep(1,nr)
  #rw<-rep(1,length(rx[,1])) #dimension change
  #nr<-length(rx[,1]) #dimension change
  rN<-length(ry)
  
  mx<-sx[sdelta==0,] #dimension change
  nm<-length(which(sdelta==0))
  mw<-rep(1,nm)
  #mw<-rep(1,length(mx[,1])) #dimension change
  #nm<-length(mx[,1]) #dimension change
  
  ###compute estimator of alpha:
  
  fn_alpha <- function(alpha){
    res <- sum((sdelta-1/(1+exp(-alpha[1]-as.matrix(sx)%*%alpha[-1]))))/n
    for (i in 1:length(sx[1,])){
      res <- c(res,sum((sdelta-1/(1+exp(-alpha[1]-as.matrix(sx)%*%alpha[-1])))*sx[,i])/n)
    }
    return(res)
  }
  
  ealpha <- nleqslv(rep(0,length(sx[1,])+1), fn_alpha)$x
  
  rpi <- 1/(1+exp(-ealpha[1]-as.matrix(rx)%*%ealpha[-1]))
  
  #mean
  etheta1 <- sum(rw*ry/rpi)/sum(rw/rpi)
  
  #median
  kcdf <- function(x) {
    temp.x <- 0.75*x-0.25*x^3+0.5
    temp.x[x>=1] <- 1
    temp.x[x<=-1] <- 0
    return(temp.x)
  }
  
  #h <- h.amise(ry, deriv.order=0, kernel="epanechnikov")[[6]]
  h <- 1.06*sd(ry)/n^(1/5)
  fq<-function(tt){
    res <- (sum(rw/rpi*kcdf((ry-tt)/h))/sum(rw/rpi)-0.5)^{2}
    return(res)
  }
  etheta2<-optimize(fq,quantile(y,c(0.05,0.95)))$minimum
  # if (abs(optimize(fq,c(min(y),max(y)))$minimum-max(y))<0.1){
  #   etheta2<-optimize(fq,c(min(y),theta0[2]+10))$minimum
  # } else {
  #   etheta2<-optimize(fq,c(min(y),max(y)))$minimum
  # }
  
  return(c(etheta1,etheta2))
}

###6. Doubly robust estimator:
fDR<-function(dat,dat.ps){
  num.col <- length(dat[1,])
  y<-dat[,1]
  x<-dat[,2:(num.col-2)] #dimension change
  delta<-dat[,num.col-1]
  sI<-dat[,num.col]
  n<-sum(sI)
  
  sx<-x[sI==1,] #dimension change
  sy<-y[sI==1]
  sw<-rep(1,length(sx[,1]))
  sdelta<-delta[sI==1]
  sN<-length(sy)
  
  rx<-sx[sdelta==1,] #dimension change
  ry<-sy[sdelta==1]
  nr<-length(which(sdelta==1))
  rw<-rep(1,nr)
  #rw<-rep(1,length(rx[,1])) #dimension change
  #nr<-length(rx[,1]) #dimension change
  rN<-length(ry)
  
  mx<-sx[sdelta==0,] #dimension change
  nm<-length(which(sdelta==0))
  mw<-rep(1,nm)
  #mw<-rep(1,length(mx[,1])) #dimension change
  #nm<-length(mx[,1]) #dimension change
  
  #PS data;
  num.col.ps <- length(dat.ps[1,])
  y.ps<-dat.ps[,1]
  x.ps<-dat.ps[,2:(num.col.ps-2)]
  delta.ps<-dat.ps[,num.col.ps-1]
  sI.ps<-dat.ps[,num.col.ps]
  n.ps<-sum(sI.ps)
  
  sx.ps<-x.ps[sI.ps==1,]
  sy.ps<-y.ps[sI.ps==1]
  sw.ps<-rep(1,length(sx.ps[,1]))
  sdelta.ps<-delta.ps[sI.ps==1]
  sN.ps<-length(sy.ps)
  
  rx.ps<-sx.ps[sdelta.ps==1,]
  ry.ps<-sy.ps[sdelta.ps==1]
  nr.ps<-length(which(sdelta.ps==1))
  rw.ps<-rep(1,nr.ps)
  #rw.ps<-rep(1,length(rx.ps[,1]))
  #nr.ps<-length(rx.ps[,1])
  rN.ps<-length(ry.ps)
  
  ###compute weighted least square estimator of beta:
  glm.dat <- data.frame(cbind(ry,rx))
  names(glm.dat) <- c("ry","X1","X2")
  model <- glm(formula = ry ~ ., data=glm.dat)
  em<-predict(model,data.frame(sx),type="response") #dimension change
  
  mem<-predict(model,data.frame(matrix(mx,ncol = 2)),type="response") #dimension change
  rem<-predict(model,data.frame(matrix(rx,ncol = 2)),type="response") #dimension change
  
  ###compute estimator of alpha:
  
  fn_alpha <- function(alpha){
    res <- sum((sdelta.ps-1/(1+exp(-alpha[1]-as.matrix(sx.ps)%*%alpha[-1]))))/n
    for (i in 1:length(sx.ps[1,])){
      res <- c(res,sum((sdelta.ps-1/(1+exp(-alpha[1]-as.matrix(sx.ps)%*%alpha[-1])))*sx.ps[,i])/n)
    }
    return(res)
  }
  
  ealpha <- nleqslv(rep(0,length(sx.ps[1,])+1), fn_alpha)$x
  
  rpi <- 1/(1+exp(-ealpha[1]-as.matrix(rx.ps)%*%ealpha[-1]))
  spi <- 1/(1+exp(-ealpha[1]-as.matrix(sx.ps)%*%ealpha[-1]))
  spi[which(spi==0)] <- 10^(-100)
  
  ###compute EL weights:
  
  re<-ry-rem
  u1<-rw*re
  elambda<-Lfv(t(u1))
  reg<-1/nr*1/(1+elambda*u1)
  reg2<-reg*rw/sum(reg*rw)
  
  ###compute the point estimate:
  
  ##matrix for eyij: (S * Sr)
  a1<-rep(em,nr)
  A1<-matrix(a1,n,nr)
  b1<-rep(re,n)
  B1<-matrix(b1,n,nr,byrow=T)
  EYIJ<-A1+B1
  wg<-sw%o%reg2
  
  #mean
  etheta1<-(sum(rw*ry/rpi)+sum((1-sdelta/spi)*rowSums(EYIJ*wg)))/sN
  
  #median
  # kcdf <- function(x) {
  #   temp.x <- as.numeric((0.75*x-0.25*x^3)>=0)
  #   temp.x[x>=1] <- 1
  #   temp.x[x<=-1] <- 0
  #   return(temp.x)
  # }
  
  kcdf <- function(x) {
    temp.x <- 0.75*x-0.25*x^3+0.5
    temp.x[x>=1] <- 1
    temp.x[x<=-1] <- 0
    return(temp.x)
  }
  
  #h <- h.amise(ry, deriv.order=0, kernel="epanechnikov")[[6]]
  h <- 1.06*sd(ry)/n^(1/5)
  fq<-function(tt){
    IEYIJ3<-matrix(kcdf((EYIJ-tt)/h),sN,nr)
    res <- ((sum(rw/rpi*kcdf((ry-tt)/h))+sum((1-sdelta/spi)*rowSums(IEYIJ3*wg)))/sN-0.5)^{2}
    return(res)
  }
  etheta2<-optimize(fq,quantile(y,c(0.05,0.95)))$minimum
  # if (abs(optimize(fq,c(min(y),max(y)))$minimum-max(y))<0.1){
  #   etheta2<-optimize(fq,c(min(y),theta0[2]+10))$minimum
  # } else {
  #   etheta2<-optimize(fq,c(min(y),max(y)))$minimum
  # }
  
  return(c(etheta1,etheta2))
}

### Jackknife;
J_CI = function(V, n, etheta0){
  etheta <- sum(V)/n
  stheta <- sum((V-etheta0)^2)/(n-1)
  
  CI_L <- etheta + qnorm(0.025)*sqrt(stheta/n)
  CI_U <- etheta + qnorm(0.975)*sqrt(stheta/n)
  # if (n>=50){
  #   CI_L <- etheta + qnorm(0.025)*sqrt(stheta/n)
  #   CI_U <- etheta + qnorm(0.975)*sqrt(stheta/n)
  # } else if(n<50){
  #   CI_L <- etheta + qt(0.025,n-1)*sqrt(stheta/n)
  #   CI_U <- etheta + qt(0.975,n-1)*sqrt(stheta/n)
  # }
  return(c(etheta,CI_L,CI_U))
}

###7. Jackknife EL
JEL_CI = function(V,n,alpha = 0.05){
  R = c()
  theta = c()
  for(i in 1:((max(V,na.rm=T)-min(V,na.rm=T))*1000)){
    theta[i] = i*.001 + min(V,na.rm=T)
    cons = function(lambda){sum((V-theta[i])/(1+lambda*(V-theta[i])))}
    d = (1-n)/n/(V-theta[i])
    upper = min(d[d>0])
    lower = max(d[d<0])
    if(cons(upper)*cons(lower)<0){
      lambda = uniroot(cons,interval = c(lower,upper))$root
    }else{
      lambda = -999
    }
    R[i] = 2*sum(log(1+lambda*(V-theta[i])))
  }
  if (n>=50){
    idx = which(R<qnorm(0.975)^2)
  } else if (n<50){
    idx = which(R<qt(0.975,n-1)^2)
  }
  #idx = which(R<1.96^2)
  return(c(theta[min(idx)],theta[max(idx)]))
}

#adjusted JEL
AJEL_CI = function(V,n,alpha = 0.05){
  R = c()
  theta = c()
  for(i in 1:((max(V,na.rm=T)-min(V,na.rm=T))*1000)){
    theta[i] = i*.001 + min(V,na.rm=T)
    V[n+1] = - log(n)/2 * mean(V[1:n],na.rm=T)
    cons = function(lambda){sum((V-theta[i])/(1+lambda*(V-theta[i])))}
    d = (1-n)/n/(V-theta[i])
    upper = min(d[d>0])
    lower = max(d[d<0])
    if(cons(upper)*cons(lower)<0){
      lambda = uniroot(cons,interval = c(lower,upper))$root
    }else{
      lambda = -999
    }
    R[i] = 2*sum(log(1+lambda*(V-theta[i])))
  }
  if (n>=50){
    idx = which(R<qnorm(0.975)^2)
  } else if (n<50){
    idx = which(R<qt(0.975,n-1)^2)
  }
  #idx = which(R<1.96^2)
  return(c(theta[min(idx)],theta[max(idx)]))
}

###8. Bayesian Jackknife EL;
BJEL_CI = function(V,n,pmu=0.5,psd=1,niter = 1e4, alpha = 0.05){
  Theta = c()
  Weight = c()
  Vmax = max(V)
  Vmin = min(V)
  for(iter in 1:niter){
    # generate theta
    theta = runif(1,Vmin,Vmax)
    cons = function(lambda){sum((V-theta)/(1+lambda*(V-theta)))}
    d = (1-n)/n/(V-theta)
    upper = min(d[d>0])
    lower = max(d[d<0])
    if(cons(upper)*cons(lower)<0){
      lambda = uniroot(cons,interval = c(lower,upper))$root
    }else{
      lambda = -999
    }
    ## calculate weight
    weight = 0
    if(lambda!=-999) weight = 1/prod(lambda*(V-theta)+1)*dnorm(theta,mean=pmu,sd = psd)
    Theta[iter] = theta
    Weight[iter] = weight
  }
  Weight = Weight/sum(Weight)
  CI = cbind(Theta,Weight,Weight)
  CI = CI[order(Theta),]
  colnames(CI) = c("theta","weight","F")
  CI[,3] = cumsum(CI[, 2])
  #theta_low = CI[max(which(CI[,3]<(alpha/2)),1),1]
  #theta_up = CI[min(which(CI[,3]>(1-alpha/2)),niter),1]
  ## find the region
  thup = 1/16
  thlow = 0
  th = thup
  id = which(CI[,2]<th)
  est = sum(CI[id,2])
  count = 0
  while(count < 100 & abs(est-alpha)>1e-3){
    count = count + 1
    if(est>alpha){
      thup = th
      th = (thup+thlow)/2
      id = which(CI[,2]<th)
      est = sum(CI[id,2])
    }else{
      thlow = th
      th = (thup+thlow)/2
      id = which(CI[,2]<th)
      est = sum(CI[id,2])
    }
  }
  id = which(CI[,2]<th)
  gap = which(diff(id)>1)
  theta_low = CI[id[gap],1]
  theta_up = CI[id[gap+1],1]
  if(length(gap)==0){
    id = setdiff(1:niter,id)
    theta_low = CI[min(id),1]
    theta_up = CI[max(id),1]
  }
  theta.id = which.max(CI[,2])
  theta = CI[theta.id,1]
  
  #return(CI)
  
  return(c(theta,theta_low,theta_up))
}

#for median;
BJEL_CI2 = function(V,n,pmu=0.5,psd=1,niter = 1e4, alpha = 0.05){
  Theta = c()
  Weight = c()
  Vmax = max(V)
  Vmin = min(V)
  for(iter in 1:niter){
    # generate theta
    theta = runif(1,Vmin,Vmax)
    cons = function(lambda){sum((V-theta)/(1+lambda*(V-theta)))}
    d = (1-n)/n/(V-theta)
    upper = min(d[d>0])
    lower = max(d[d<0])
    if(cons(upper)*cons(lower)<0){
      lambda = uniroot(cons,interval = c(lower,upper))$root
    }else{
      lambda = -999
    }
    ## calculate weight
    weight = 0
    if(lambda!=-999) weight = 1/prod(lambda*(V-theta)+1)*dnorm(theta,mean=pmu,sd = psd)
    Theta[iter] = theta
    Weight[iter] = weight
  }
  Weight = Weight/sum(Weight)
  CI = cbind(Theta,Weight,Weight)
  CI = CI[order(Theta),]
  colnames(CI) = c("theta","weight","F")
  CI[,3] = cumsum(CI[, 2])
  #theta_low = CI[max(which(CI[,3]<(alpha/2)),1),1]
  #theta_up = CI[min(which(CI[,3]>(1-alpha/2)),niter),1]
  ## find the region
  thup = 1/16
  thlow = 0
  th = thup
  id = which(CI[,2]<th)
  est = sum(CI[id,2])
  count = 0
  while(count < 100 & abs(est-alpha)>1e-3){
    count = count + 1
    if(est>alpha){
      thup = th
      th = (thup+thlow)/2
      id = which(CI[,2]<th)
      est = sum(CI[id,2])
    }else{
      thlow = th
      th = (thup+thlow)/2
      id = which(CI[,2]<th)
      est = sum(CI[id,2])
    }
  }
  id = which(CI[,2]<th)
  gap = which(diff(id)>1)
  theta_low = CI[id[gap],1]
  theta_up = CI[id[gap+1],1]
  if(length(gap)==0){
    id = setdiff(1:niter,id)
    theta_low = CI[min(id),1]
    theta_up = CI[max(id),1]
  }
  # theta.id = which.max(CI[,2])
  # theta = CI[theta.id,1]
  theta=weighted.median(CI[,1],CI[,2])
  
  #return(CI)
  
  return(c(theta,theta_low,theta_up))
}

###9. Adjusted Bayesian Jackknife EL;
ABJEL_CI = function(V,n,pmu=0.5,psd=1,niter = 1e4,alpha = 0.05){
  Theta = c()
  Weight = c()
  Vmax = max(V)
  Vmin = min(V)
  for(iter in 1:niter){
    # generate theta
    #theta = rnorm(1,mean = .8183, sd = 1)
    theta = runif(1,Vmin,Vmax)
    
    V[n+1] = - log(n)/2 * mean(V[1:n]) + theta + log(n)/2 * theta
    cons = function(lambda){sum((V-theta)/(1+lambda*(V-theta)))}
    d = (1-n)/n/(V-theta)
    upper = min(d[d>0])
    lower = max(d[d<0])
    if(cons(upper)*cons(lower)<0){
      lambda = uniroot(cons,interval = c(lower,upper))$root
    }else{
      lambda = -999
    }
    ## calculate weight
    weight = 0
    #if(lambda!=-999) weight = 1/prod(lambda*(V-theta)+1)*dnorm(theta,mean = pmu, sd = psd)/dnorm(theta,mean = .8183, sd = 1)
    if(lambda!=-999) weight = 1/prod(lambda*(V-theta)+1)*dnorm(theta,mean = pmu, sd = psd)
    Theta[iter] = theta
    Weight[iter] = weight
  }
  Weight = Weight/sum(Weight)
  CI = cbind(Theta,Weight,Weight)
  CI = CI[order(Theta),]
  colnames(CI) = c("theta","weight","F")
  CI[,3] = cumsum(CI[, 2])
  #theta_low = CI[max(which(CI[,3]<(alpha/2)),1),1]
  #theta_up = CI[min(which(CI[,3]>(1-alpha/2)),niter),1]
  ## find the region
  thup = 1/16
  thlow = 0
  th = thup
  id = which(CI[,2]<th)
  est = sum(CI[id,2])
  count = 0
  while(count < 100 & abs(est-alpha)>1e-3){
    count = count + 1
    if(est>alpha){
      thup = th
      th = (thup+thlow)/2
      id = which(CI[,2]<th)
      est = sum(CI[id,2])
    }else{
      thlow = th
      th = (thup+thlow)/2
      id = which(CI[,2]<th)
      est = sum(CI[id,2])
    }
  }
  id = which(CI[,2]<th)
  gap = which(diff(id)>1)
  theta_low = CI[id[gap],1]
  theta_up = CI[id[gap+1],1]
  if(length(gap)==0){
    id = setdiff(1:niter,id)
    theta_low = CI[min(id),1]
    theta_up = CI[max(id),1]
  }
  theta.id = which.max(CI[,2])
  theta = CI[theta.id,1]
  
  return(c(theta,theta_low,theta_up))
}

ABJEL_CI2 = function(V,n,pmu=0.5,psd=1,niter = 1e4,alpha = 0.05){
  Theta = c()
  Weight = c()
  Vmax = max(V)
  Vmin = min(V)
  for(iter in 1:niter){
    # generate theta
    #theta = rnorm(1,mean = .8183, sd = 1)
    theta = runif(1,Vmin,Vmax)
    
    V[n+1] = - log(n)/2 * mean(V[1:n]) + theta + log(n)/2 * theta
    cons = function(lambda){sum((V-theta)/(1+lambda*(V-theta)))}
    d = (1-n)/n/(V-theta)
    upper = min(d[d>0])
    lower = max(d[d<0])
    if(cons(upper)*cons(lower)<0){
      lambda = uniroot(cons,interval = c(lower,upper))$root
    }else{
      lambda = -999
    }
    ## calculate weight
    weight = 0
    #if(lambda!=-999) weight = 1/prod(lambda*(V-theta)+1)*dnorm(theta,mean = pmu, sd = psd)/dnorm(theta,mean = .8183, sd = 1)
    if(lambda!=-999) weight = 1/prod(lambda*(V-theta)+1)*dnorm(theta,mean = pmu, sd = psd)
    Theta[iter] = theta
    Weight[iter] = weight
  }
  Weight = Weight/sum(Weight)
  CI = cbind(Theta,Weight,Weight)
  CI = CI[order(Theta),]
  colnames(CI) = c("theta","weight","F")
  CI[,3] = cumsum(CI[, 2])
  #theta_low = CI[max(which(CI[,3]<(alpha/2)),1),1]
  #theta_up = CI[min(which(CI[,3]>(1-alpha/2)),niter),1]
  ## find the region
  thup = 1/16
  thlow = 0
  th = thup
  id = which(CI[,2]<th)
  est = sum(CI[id,2])
  count = 0
  while(count < 100 & abs(est-alpha)>1e-3){
    count = count + 1
    if(est>alpha){
      thup = th
      th = (thup+thlow)/2
      id = which(CI[,2]<th)
      est = sum(CI[id,2])
    }else{
      thlow = th
      th = (thup+thlow)/2
      id = which(CI[,2]<th)
      est = sum(CI[id,2])
    }
  }
  id = which(CI[,2]<th)
  gap = which(diff(id)>1)
  theta_low = CI[id[gap],1]
  theta_up = CI[id[gap+1],1]
  if(length(gap)==0){
    id = setdiff(1:niter,id)
    theta_low = CI[min(id),1]
    theta_up = CI[max(id),1]
  }
  # theta.id = which.max(CI[,2])
  # theta = CI[theta.id,1]
  theta=weighted.median(CI[,1],CI[,2])
  
  return(c(theta,theta_low,theta_up))
}
