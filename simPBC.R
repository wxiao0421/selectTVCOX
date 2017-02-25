##########################################################################################
# Numerical study with PBC data 
##########################################################################################
rm(list=ls())
library(survival)
library(penalized)
source("KGNG.R")
set.seed(919889)

##########################################################################################
# Load in the data 
##########################################################################################
# read in the data (17 covariates) 
x <- as.matrix(read.table("xliver.dat"))
dimnames(x)[[2]] <- names(pbc)[-1:-3]
y <- as.matrix(read.table("yliver.dat"))
dimnames(y)[[2]]<-c("days", "status")
data<-list(CenTime=as.vector(y[,1]),Delta=as.vector(y[,2]),Zmat=x)
data$Delta<-as.numeric(data$Delta==2)
# take log on serum bilirubin, albumin and prothrombin
data$Zmat[,c(8,10,16)]<-log(data$Zmat[,c(8,10,16)])
# rescale the unites of age to 10 years
data$Zmat[,2]<-data$Zmat[,2]/10
# rescale copper to mg/day
data$Zmat[,11]<-data$Zmat[,11]/1000
# randomize data
rands<-runif(length(data$CenTime),0,1)
orders<-order(rands)
data$Zmat<-data$Zmat[orders,]
data$CenTime<-data$CenTime[orders]
data$Delta<-data$Delta[orders]
# standardize data to data2
data2<-data
data2$Zmat<-scale(data$Zmat, center=TRUE, scale=TRUE)
a1<-attr(data2$Zmat,"scaled:center")
b1<-attr(data2$Zmat,"scaled:scale")

# set up parameters
n<-length(data2$Delta)
rbound=3200
hvect<-2000 # fix bandwidth to 2000
out.variance<-T
kappa<-NA
thetaVect<-exp(seq(-6,-4,0.5))*n
#(theta1,theta2)=thetaMatrix=n*exp(thetaMatrix1)
theta1<-seq(-6,-4,0.5)
theta2<-seq(-6,-4,0.5)
thetaMatrix1<-matrix(NA,length(theta1)*length(theta2),2)
k<-0
for(i in 1:length(theta1))
{
  for(j in 1:length(theta2))
  {
    k<-k+1
    thetaMatrix1[k,]<-c(theta1[i],theta2[j])
  }
}
thetaMatrix<-exp(thetaMatrix1)*n

# fit KGNG
prelim.step <- F
Nmin=0
KGNG.out <- KGNG(data2, hvect, prelim.step, out.variance, kappa, rbound, 
                 thetaVect=thetaVect, thetaMatrix=thetaMatrix)

# fit KGNG2
prelim.step <- T
KGNG2.out <- KGNG(data2, hvect, prelim.step, out.variance, kappa, rbound,
                  thetaVect=thetaVect, thetaMatrix=thetaMatrix)

save(data, data2, a1, b1, KGNG.out, KGNG2.out, file="Result.RData")

##########################################################################################
# Rescale the result back to the original scale and print the result 
##########################################################################################
load("Result.RData")

# coefficient functions (KGNG) 
print(diag(1/b1)%*%KGNG.out$coefMatrix)
# variances of time-varying effect coefficient functions (KGNG)   
print(sweep(KGNG.out$SE_NC,1,b1,FUN="/"))
# variances of constant effect coefficients (KGNG) 
print(KGNG.out$SE_C/b1[KGNG.out$select==1])

# coefficient functions (KGNG2) 
print(diag(1/b1)%*%KGNG2.out$coefMatrix)
# variances of time-varying effect coefficient functions (KGNG2)   
print(sweep(KGNG2.out$SE_NC,1,b1,FUN="/"))
# variances of constant effect coefficients (KGNG2) 
print(KGNG2.out$SE_C/b1[KGNG2.out$select==1])










