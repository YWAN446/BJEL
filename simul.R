#packages
library(Hmisc)
library(MASS)
library(kedd)
library(nleqslv)
library(spatstat)
library(boot)
###Simulation###
#scenarios

#load functions
source("./func.R")
#parameters
seed<-c(15,14,13,12,11)
B=2000
N=10000
n.sample=50
alpha=0.05
niter=10000

sigma1 <- 1
sigma2 <- sqrt(n.sample)

#generate data;
dat1_100<-gedata(B=B,N=N,n=n.sample,seed=seed)
Dat0 = dat1_100[[1]]
theta0 = dat1_100[[2]]
rm(dat1_100)

#correctly specified model as illustration;
ind1 <- 4
ind2 <- 4
Dat <- Dat0[,-ind1,]

#only run for iteration 1 for illustration;
i=1
#Use SFI estimator as an example; Can replace fSFIS with fPS for PS estimator.
tn = fSFIS(Dat[,,i])
n = length(which(Dat[,5,i]==1))
V1 <- c()
V2 <- c()
for (j in 1:n){
  tnj <- fSFIS(Dat[-which(Dat[,5,i]==1)[j],,i])
  V1[j] = (n*tn - (n-1)*tnj)[1]
  V2[j] = (n*tn - (n-1)*tnj)[2]
}
#mean
Res[i,1] <- tn[1]
Res[i,2:4] <- J_CI(V1,n,tn[1])
Res[i,5:6] <- JEL_CI(V1,n,alpha = alpha)
Res[i,7:8] <- AJEL_CI(V1,n,alpha = alpha)
Res[i,9:11] <- BJEL_CI(V1,n,pmu=theta0[1],psd=sigma1,niter = niter, alpha = alpha)
Res[i,12:14] <- ABJEL_CI(V1,n,pmu=theta0[1],psd=sigma1,niter = niter, alpha = alpha)
Res[i,15:17] <- BJEL_CI(V1,n,pmu=theta0[1],psd=sigma2,niter = niter, alpha = alpha)
Res[i,18:20] <- ABJEL_CI(V1,n,pmu=theta0[1],psd=sigma2,niter = niter, alpha = alpha)
#median
Res[i,21] <- tn[2]
Res[i,22:24] <- J_CI(V2,n,tn[2])
Res[i,25:26] <- JEL_CI(V2,n,alpha = alpha)
Res[i,27:28] <- AJEL_CI(V2,n,alpha = alpha)
Res[i,29:31] <- BJEL_CI2(V2,n,pmu=theta0[2],psd=sigma1,niter = niter, alpha = alpha)
Res[i,32:34] <- ABJEL_CI2(V2,n,pmu=theta0[2],psd=sigma1,niter = niter, alpha = alpha)
Res[i,35:37] <- BJEL_CI2(V2,n,pmu=theta0[2],psd=sigma2,niter = niter, alpha = alpha)
Res[i,38:40] <- ABJEL_CI2(V2,n,pmu=theta0[2],psd=sigma2,niter = niter, alpha = alpha)

#For DR estimator;
Dat1 <- Dat0[,-ind1,]
Dat2 <- Dat0[,-ind2,]

tn = fDR(Dat1[,,i],Dat2[,,i])
n = length(which(Dat1[,5,i]==1))
V1 <- c()
V2 <- c()
for (j in 1:n){
  tnj <- fDR(Dat1[-which(Dat1[,5,i]==1)[j],,i],Dat2[-which(Dat2[,5,i]==1)[j],,i])
  V1[j] = (n*tn - (n-1)*tnj)[1]
  V2[j] = (n*tn - (n-1)*tnj)[2]
}
#mean
Res[i,1] <- tn[1]
Res[i,2:4] <- J_CI(V1,n,tn[1])
Res[i,5:6] <- JEL_CI(V1,n,alpha = alpha)
Res[i,7:8] <- AJEL_CI(V1,n,alpha = alpha)
Res[i,9:11] <- BJEL_CI(V1,n,pmu=theta0[1],psd=sigma1,niter = niter, alpha = alpha)
Res[i,12:14] <- ABJEL_CI(V1,n,pmu=theta0[1],psd=sigma1,niter = niter, alpha = alpha)
Res[i,15:17] <- BJEL_CI(V1,n,pmu=theta0[1],psd=sigma2,niter = niter, alpha = alpha)
Res[i,18:20] <- ABJEL_CI(V1,n,pmu=theta0[1],psd=sigma2,niter = niter, alpha = alpha)
#median
Res[i,21] <- tn[2]
Res[i,22:24] <- J_CI(V2,n,tn[2])
Res[i,25:26] <- JEL_CI(V2,n,alpha = alpha)
Res[i,27:28] <- AJEL_CI(V2,n,alpha = alpha)
Res[i,29:31] <- BJEL_CI2(V2,n,pmu=theta0[2],psd=sigma1,niter = niter, alpha = alpha)
Res[i,32:34] <- ABJEL_CI2(V2,n,pmu=theta0[2],psd=sigma1,niter = niter, alpha = alpha)
Res[i,35:37] <- BJEL_CI2(V2,n,pmu=theta0[2],psd=sigma2,niter = niter, alpha = alpha)
Res[i,38:40] <- ABJEL_CI2(V2,n,pmu=theta0[2],psd=sigma2,niter = niter, alpha = alpha)
