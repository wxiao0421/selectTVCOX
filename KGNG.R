library(survival)
library(penalized)

##
## @fn KGNG(data, hvect, prelim.step, out.variance, kappa, rbound, nlength, thetaVect, thetaMatrix, positive, Nmin, Mpermu)
## @brief calculate KGNG (kernel group nonnegative garrote estimator) or KGNG2 (KGNG with 
##        preliminary step) in time-varying coefficient Cox model. This method can automatically  
##        identify the structure of covariate, i.e., covariates with null effect(O),  
##        constant effect (C) and truly time-varying effect (NC), and estimate the corresponding  
##        regression coefficients.
## @param data, a list of three components, Zmat (covariates, n x p matrix), 
##        Delta (censoring indicator, n x 1 vector, where Delta[i]=1(Ti<=Ci)) and 
##        CenTime (observed survival times, n x 1 vector, where CenTime[i]=min(Ti, Ci)). 
##        Here Ti and Ci are the failure time and the censored time of the ith subject, respectively.
## @param hvect vector of bandwidth, length n. When n>1, a K-fold CV is automatically called 
##        to find the optimal bandwidth. When n=1, the optimal bandwidth is hvect.
## @param prelim.step logical value: if TRUE, KGNG2 is applied. Otherwise, KGNG is applied.
## @param out.variance logical value: if TRUE, variances are calculated.
## @param kappa integer, the number of split parts in K-fold cross validation.
## @param rbound numeric, tau value.
## @param nlength integer, the number of grid points on time line.
## @param thetaVect vector, the grid points of tuning parameter in Step 0.
## @param thetMatrix matrix, dim K x 2, the grid points of tuning parameters in Step 3.
## @param positive logical vaue: if TRUE, lambda is restricted to be positive.
## @param Nmin integer, minimal observed failure times in a window when maximize the local 
##        partial likelihood.
## @param Mpermu integer, number of MC samples to calculate the cutoff point of confidence band
## @return list with the following key-->value pairs:
##         select-->p x 1 vector, structure of the covariates. 0 represents covariate with 
##                  null effect, 1 represents covariate with constant effect and 2 represents   
##                  covariate with time-varying effect.
##         select.nonzero-->p x 1 vector, structure of the covariates. 0 represents covariate  
##                  with null effect, and 1 represents covariate with non-null effect.
##         coefMatrix.ini-->initial estimats of coefficient matrix with dimension p x M, 
##                          where M=nlength.
##         coefMatrix-->final estimats of coefficient matrix with dimension p x M, where 
##                      M=nlength.
##         theta-->3 x 1 vector. theta[1] is the tuning parameter of Step 0; (theta[2], theta[3])
##         are the tuning parameters of Step 3.
##         SE_C-->standard errors of covariates with constant effect
##         cov_C-->covariance matrix of covariates with constant effect
##         SE_NC-->p x M matrix, standard errors of covariates with time-varying effect
##         c_alpha-->p x 1 vector, cutoff point of confidence band 
##         time.h-->running time to find the optimal bandwidth with K-fold CV 
##         time.prelim-->running time for Step 0 
##         time.final-->running time for Step 1 and 2 
##
KGNG <- function(data, hvect, prelim.step, out.variance, kappa=5, rbound=NA, nlength=100, 
                 thetaVect=NA, thetaMatrix=NA, positive=TRUE, Nmin=0, Mpermu=10000){
  
  #set up all parameters
  n = nrow(data$Zmat)
  p = ncol(data$Zmat)
  lbound = 0
  if (length(thetaVect)==0){
    if (is.na(thetaVect)) thetaVect = exp(seq(0,-8,-0.5))*n
  }
  if (length(thetaMatrix)==0){
    if (is.na(thetaMatrix)) {
      theta1<-seq(-8,0,1)
      theta2<-seq(-8,0,1)
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
    }
  }
  if (is.na(rbound)) rbound = max(data$CenTime) else {
    # truncate data set by rbound
    data.copy<-data
    data$Delta[data.copy$CenTime>rbound]<-0
    data$CenTime[data.copy$CenTime>rbound]<-rbound 
  }
  timePoint = seq(lbound, rbound, length.out=nlength) 
  # calculate the initial values of beta by fitting PH model
  phfit <- coxph(Surv(CenTime, Delta)~Zmat,data)
  beta0.ini <- phfit$coef
  coefini =  beta0.ini%*%t(rep(1, nlength))
  
  # find optimal bandwidth
  ptm <- proc.time()
  if(length(hvect)==1) h = hvect else {
    optH.out <- searchBandwidth(kappa, data, hvect, estBetaTime, 1:p, beta0.ini)
    h <- hvect[which.min(optH.out$APE)]
  }
  # Stop the clock
  time.h = (proc.time() - ptm)[3]
  ptm <- proc.time()  
  
  if (prelim.step==TRUE){
    prelim.out <- prelimStep(data, n, p, beta0.ini, nlength, timePoint, Nmin, h, 
                             lbound, rbound, thetaVect, positive)
    Index.nonzero = prelim.out$Index.nonzero
    theta.star.opt = prelim.out$theta.star.opt
    time.prelim = prelim.out$time.prelim
    coefMatrix = prelim.out$coefMatrix
    coefMatrix2 = prelim.out$coefMatrix2
  } else {
    Index.nonzero = 1:p
    theta.star.opt = NA
    time.prelim = NA
    coefMatrix = NA
    coefMatrix2 = NA    
  }  
  select.nonzero <- rep(0,p)
  if (length(Index.nonzero)!=0){
    select.nonzero[Index.nonzero] <- 1
  }
  # get data2 where only important covariates are kept 
  # or data2=data if prelim.Step=FALSE
  data2 <- data
  data2$Zmat <- as.matrix(data$Zmat[,Index.nonzero])
  # Stop the clock
  time.prelim = (proc.time() - ptm)[3]
  ptm <- proc.time() 
  
  ### selection between zero, constant and non-constant covariates ### 
  ### the group nonnegative garrote step ###
  if (length(Index.nonzero)!=0){
    # use the same bandwidth for step 2 and step 1
    h2=h    
    # estimate coefficients with local constant partial likelihood
    if (length(Index.nonzero)>1) coefini.Index <- coefini[Index.nonzero,] else 
      coefini.Index <- as.matrix(t(coefini[Index.nonzero,]))
    const.fit <- estBeta_h(h2, data2, timePoint, Nmin, coefini.Index)
    coefMatrix.const <- const.fit$betaMatrix
    # betaV2=Sigma(t)
    betaV2 <- const.fit$betaV
    coefMatrix.const1 <- matrix(0, p, nlength)
    coefMatrix.const1[Index.nonzero,] <- coefMatrix.const
    # estimate beta mean for second step
    p.const <- length(Index.nonzero)
    betaMean.const <- rowMeans(coefMatrix.const)
    star.const <- sweep(coefMatrix.const,1,betaMean.const,FUN="-")
    
    if(out.variance==TRUE){
      # calculate the variance covariance matrix of Beta tilde and m tilde, and Beta tilde's 
      # corresponding 95% CI and CB
      ini.KGNG2.object <-
        iniCI(const.fit, data2, coefMatrix.const, timePoint, rbound, lbound, n, p.const, h2, 
              Mpermu)
      invbetaV2 <- ini.KGNG2.object$invbetaV
      # variance covariance matrix of m tilde
      sigmaM2 <- ini.KGNG2.object$sigmaM
      # SE of m tilde
      SEM2 <- ini.KGNG2.object$SEM
      # SE of tilde beta(t)
      SEBeta2 <- ini.KGNG2.object$SEBeta
      # constant of 95% confidence band
      c_alpha2 <- ini.KGNG2.object$c_alpha
    }
    
    
    # get Y and X matrix of Cox model with time-dependent covariates for step three
    k <- 0
    for(i in 1:n)
    {
      for (j in 1:(nlength-1))
      {
        if (timePoint[j]<data2$CenTime[i]) k <- k+1  
      }
    }   
    start <- rep(NA, k)
    stop <- rep(NA, k)
    cen.time<-rep(NA, k)
    X <- matrix(NA, k, (p.const*2))  
    k2 <- 0
    for(i in 1:n)
    {
      for (j in 1:(nlength-1))
      {
        if (timePoint[j]<data2$CenTime[i]) 
        {
          k2 <- k2+1
          start[k2] <- timePoint[j]
          stop[k2] <- timePoint[(j+1)]
          if(timePoint[(j+1)]>=data2$CenTime[i]) cen.time[k2]=data2$Delta[i] else cen.time[k2]=0
          X[k2,1:p.const] <- data2$Zmat[i,]*betaMean.const
          X[k2,(p.const+1):(2*p.const)] <- 
            data2$Zmat[i,]*((coefMatrix.const[,j]+coefMatrix.const[,(j+1)])/2-betaMean.const)
        }  
      }
    }  
    Y <- Surv(start,stop,cen.time) #Define Y as a survival subject    
    # group nonnegative garrote      
    for(i in 1:nrow(thetaMatrix))
    {
      theta_fixed <- thetaMatrix[i,]
      theta1_fixed <- theta_fixed[1]
      theta2_fixed <- theta_fixed[2]
      ratio_fixed <- theta2_fixed/theta1_fixed
      X2 <- X
      X2[,(p.const+1):(2*p.const)] <- X[,(p.const+1):(2*p.const)]/ratio_fixed
      # use "penalized" function to get lambda_hat
      pen <- penalized(Y, penalized = X2, lambda1=theta1_fixed, positive=positive) 
      lambda_fixed <- coefficients(pen,"penalized")
      lambda1_fixed <- lambda_fixed[1:p.const]
      lambda2_fixed <- lambda_fixed[(p.const+1):(2*p.const)]/ratio_fixed
      select2.2Dim.temp<-rep(0,p.const)
      for (j in 1:p.const)
      {
        if (lambda2_fixed[j]!=0) select2.2Dim.temp[j]=2 else {
          if (lambda1_fixed[j]!=0) select2.2Dim.temp[j]=1 else select2.2Dim.temp[j]=0    
        }
      }  
      select.const.temp <- rep(0,p)
      select.const.temp[Index.nonzero] = select2.2Dim.temp        
      # etimated coefficients for third step with fixed theta  
      coefMatrix.const2.temp <- matrix(0, p, nlength)
      if (length(Index.nonzero)>1)
        coefMatrix.const2.temp[Index.nonzero,] <- 
        sweep(diag(lambda2_fixed)%*%star.const, 1, betaMean.const*lambda1_fixed, FUN="+")
      else if (length(Index.nonzero)==1)
        coefMatrix.const2.temp[Index.nonzero,] <- 
        sweep(diag(lambda2_fixed,1,1)%*%star.const, 1, betaMean.const*lambda1_fixed, FUN="+")    
      #Calculate BIC value
      logLike <- loglik(pen)
      df1 <- sum(lambda1_fixed!=0)
      df2 <- sum(lambda2_fixed!=0)
      bic <- -2*logLike+df1*log(n)+df2*log(n*h2/(rbound-lbound))/(h2/(rbound-lbound))
      if (i==1){
        bic.OPT <- bic
        select.const <- select.const.temp
        coefMatrix.const2 <- coefMatrix.const2.temp
        select2.2Dim <- select2.2Dim.temp
        theta1.opt = theta1_fixed
        theta2.opt = theta2_fixed 
      } else if (bic<bic.OPT){
        bic.OPT <- bic
        select.const <- select.const.temp
        coefMatrix.const2 <- coefMatrix.const2.temp
        select2.2Dim <- select2.2Dim.temp
        theta1.opt = theta1_fixed
        theta2.opt = theta2_fixed 
      }        
    }
    
    if(out.variance==TRUE){
      #calculate SEBeta_m2
      order2 <- order(select2.2Dim)
      p1 <- sum(select2.2Dim==0)
      p2 <- sum(select2.2Dim==1)
      p3 <- sum(select2.2Dim==2)
      betaMean_order <- betaMean.const[order2]
      star_order <- star.const[order2,]
      if (p2>=1){
        m_C <- betaMean_order[(p1+1):(p1+p2)]
        if (p2>1){
          diag_C = diag(m_C)
        } else diag_C = m_C
        #standard error of hat beta_m
        SEBeta2_m <- rep(0, p2)
        #variance covariance matrix of hat beta_m
        sigmaBeta2_m <- matrix(0, p2, p2)
        if (p3>=1){
          m_NC <- betaMean_order[(p1+p2+1):p.const]
          star_NC <- star_order[(p1+p2+1):p.const,]
        }
        for (i in 1:nlength){
          temp1 <- matrix(0, p2+2*p3, p.const)
          A5 <- matrix(0, p2+2*p3, p2+2*p3)
          for(j in 1:nlength){
            P1A0s<-matrix(0, (p2+2*p3), p.const)
            P1A0s[1:p2, (p1+1):(p1+p2)] <- diag_C
            if (p3>=1){
              P1A0s[(p2+1):(p2+p3),(p1+p2+1):p.const] <- diag(m_NC)
              P1A0s[(p2+p3+1):(p2+2*p3),(p1+p2+1):p.const] <- diag(star_NC[,j])
            }
            sigma_s <- matrix(betaV2[,j], p.const, p.const)
            B3_s <- matrix(0,p.const,p2+2*p3)
            B3_s[(p1+1):(p1+p2),1:p2] <- diag_C
            if (p3>=1){
              B3_s[(p1+p2+1):p.const,(p2+1):(p2+p3)] <- diag(m_NC)
              B3_s[(p1+p2+1):p.const,(p2+p3+1):(p2+2*p3)] <- diag(star_NC[,j])
            }
            temp1 <- temp1+P1A0s%*%sigma_s/nlength
            A5 <- A5+P1A0s%*%sigma_s%*%B3_s*(rbound-lbound)/nlength
          }
          P1A0u <- matrix(0, (p2+2*p3), p.const)
          P1A0u[1:p2, (p1+1):(p1+p2)] <- diag_C
          if (p3>=1){
            P1A0u[(p2+1):(p2+p3), (p1+p2+1):p.const] <- diag(m_NC)
            P1A0u[(p2+p3+1):(p2+2*p3),(p1+p2+1):p.const] <- diag(star_NC[,i])
          }
          sigma_u <- matrix(betaV2[,i],p.const,p.const)
          invsigma_u <- matrix(invbetaV2[,i],p.const,p.const)
          B5 <- B4 <- matrix(0,p.const,p.const)
          B5[(p1+1):p.const, (p1+1):p.const] <- diag(p2+p3)
          if (p3>=1){
            B4[(p1+p2+1):p.const, (p1+p2+1):p.const] <- diag(p3)
          }
          temp2 <- P1A0u-temp1%*%B5%*%invsigma_u-P1A0u%*%sigma_u%*%B4%*%invsigma_u
          B6 <- matrix(0,p2,p2+2*p3)
          B6[1:p2,1:p2] <- diag(p2)
          D_u <- diag_C%*%B6%*%solve(A5,temp2)
          B7 <- matrix(0,p2,p.const)
          B7[1:p2, (p1+1):(p1+p2)] <- diag(p2)
          D_addition <- 1/(rbound-lbound)*B7%*%invsigma_u
          V_u <- (D_u+D_addition)%*%sigma_u%*%t(D_u+D_addition)
          sigmaBeta2_m <- sigmaBeta2_m+V_u*(rbound-lbound)/nlength/n
        }
        SEBeta2_m <- sqrt(diag(sigmaBeta2_m))
        #         print(SEBeta2_m)
        #         print(SEM2[select2.2Dim==1])
        #         eigen(sigmaM2[select2.2Dim==1,select2.2Dim==1]-sigmaBeta2_m)
      } else {
        SEBeta2_m = sigmaBeta2_m =NA  
      }
    } else {
      SEBeta2_m = sigmaBeta2_m = NA   
    }
  } else {
    coefMatrix.const1 = coefMatrix.const2 = coefMatrix2
    select.const = select.nonzero
    theta1.opt = theta2.opt = NA
    h2=NA
    SEBeta2_m = sigmaBeta2_m = NA
  }
  
  # extend estimates of variance of time-varying components to the full dimension
  if (out.variance==TRUE){
    SEBeta2full <- matrix(NA, p, nlength)
    SEBeta2full[Index.nonzero,] <- SEBeta2
    #SEBeta2full[select.const!=2,]<-NA
    c_alpha2full <- rep(NA, p)
    c_alpha2full[Index.nonzero] <- c_alpha2
    #c_alpha2full[select.const!=2]<-NA
  } else {
    SEBeta2full <- matrix(NA, p, nlength)
    c_alpha2full <- rep(NA, p)
  }
  
  # Stop the clock
  time.final = (proc.time() - ptm)[3]
  
  #   # plot beta(t)
  #   x11()
  #   par(mfrow=c(1,1))
  #   matplot(timePoint, t(coefMatrix.const2), type="l", xlab="t", ylab="beta(t)", 
  #           main="Estimates of Coefficient Functions with Kernel Group Nonnegative Garrote")
  
  return(list(select=select.const, select.nonzero=select.nonzero, 
              coefMatrix.ini=coefMatrix.const1, coefMatrix=coefMatrix.const2, 
              h=h, theta=c(theta.star.opt, theta1.opt, theta2.opt),
              SE_C=SEBeta2_m, cov_C=sigmaBeta2_m, SE_NC=SEBeta2full, c_alpha=c_alpha2full, 
              time.h=time.h, time.prelim=time.prelim, time.final=time.final))  
}

##
## @fn EpanKernel(x)
## @brief Epanechnikov Kernel
## @return K(x)
##
Ginv = function(x, tolerance = 1e-5)
{
  x.svd <- svd(x)
  d <-x.svd$d
  index <- d>=tolerance
  ds <- diag(1/x.svd$d[index])
  u <- x.svd$u
  v <- x.svd$v
  us <- as.matrix(u[, index])
  vs <- as.matrix(v[, index])
  x.ginv <- vs %*% ds %*% t(us)
  return(x.ginv)
}

##
## @fn EpanKernel(x)
## @brief Epanechnikov Kernel
## @return K(x)
##
EpanKernel = function(x)
{
  if(x>=-1 && x<=1) y=3*(1-x^2)/4 else y=0
  return(y)
}

##
## @fn loclikelihood1(beta, t, h, Zmat, CenTime, Delta)
## @brief calculate the local constant partial likelihood ll(beta) at time t with 
##        bandwidth h. 
## @param beta p x 1 vector, coefficient vector
## @param t numeric, time point
## @param h numeric, bandwidth
## @param Zmat n x p matrix, covariate matrix
## @param CenTime n x 1 vector, observed survival times (CenTime[i]=min(Ti, Ci))
## @param Delta n x 1 vector, censoring indicators (Delta[i]=1(Ti<=Ci))
## @return negll(beta) numeric, negative local constant partial likelihood
##
loclikelihood1 = function(beta, t, h, Zmat, CenTime, Delta)
{
  p <- ncol(Zmat)
  n <- length(Delta)
  # set initial value 
  ll <- 0
  llprime <- rep(0,p)
  # get minus loclikelihood and its gradient
  for (i in which(Delta==1 & abs((CenTime-t)/h)<1))
  {
    K <- EpanKernel((CenTime[i]-t)/h)
    M <- exp(as.vector(t(beta)%*%t(Zmat)))
    Y <- as.numeric(CenTime>=CenTime[i])
    G <- as.numeric(t(Y)%*%M)
    YM <- Y*M
    Gprime <- as.vector(t(YM)%*%Zmat)
    ll <- ll+as.numeric(K*(t(beta)%*%Zmat[i,]-log(G)))
    llprime <- llprime+K*(Zmat[i,]-Gprime/G)
  }
  negll <- -ll
  attr(negll, "gradient") <- -llprime    
  return(negll)
}

##
## @fn estBeta_h(h, data, timePoint, Nmin, coefTrue)
## @brief calculate the local constant partial maximum likelihood estimates with bandwidth 
##        h. 
## @param h numeric, bandwidth
## @param data, a list of Zmat (covariates matrix), Delta (censoring indicator, 
##        Delta[i]=1(Ti<=Ci)) and CenTime (observed survival times, CenTime[i]=min(Ti, Ci)), 
##        where Ti is the failure time of ith subject and Ci is the censored time of ith 
##        subject.
## @param timepoint kn x 1 vector, grid time points to estimate regression coefficients of 
##        time-varying coefficient Cox model
## @param Nmin integer, each local window has at least Nmin observed failure time points
## @param coefIni p x kn matrix, initial values of regression coefficients, where 
##        coefIni[,k] is the initial value of beta at timepoint t=timePoint[k] 
## @return list with the following key-->value pairs:
##         betaMatrix-->estimated coefficient matrix with dimension p x kn
##         betaV-->I(beta, t)
##         countV-->kn x 1 vector of counts of observed survival times within each windows 
##
estBeta_h = function(h, data, timePoint, Nmin, coefIni)
{
  Delta <- data$Delta
  CenTime <- data$CenTime
  Zmat <- data$Zmat
  n <- nrow(Zmat)
  p <- ncol(Zmat)
  ntimes <- length(timePoint)
  betaMatrix <- matrix(NA, p, ntimes)
  betaV <- matrix(NA, p*p, ntimes)
  countV <- rep(NA, ntimes)
  for (i in 1:ntimes)
  {
    # choose a bandwidth with at least Nmin observed failure time at each window
    h1 = h
    count <- sum(CenTime<=(timePoint[i]+h1) & CenTime>=(timePoint[i]-h1) & Delta==1)
    while(count<Nmin) 
    {
      h1 = h1*1.5;
      count <- sum(CenTime<=(timePoint[i]+h1) & CenTime>=(timePoint[i]-h1) & Delta==1)
    }     
    beta <- nlm(loclikelihood1, coefIni[,i], t=timePoint[i], 
                  h=h1, Zmat=Zmat, CenTime=CenTime,  Delta=Delta)$estimate
    betaMatrix[,i] <- beta
    countV[i] <- count
    #caluclate I(beta, t)
    I <- matrix(0, p, p)
    for (j in which(Delta==1 & abs((CenTime-timePoint[i])/h1)<1))
    {
      t <- timePoint[i]
      K <- EpanKernel((CenTime[j]-t)/h1)/h1
      M <- exp(as.vector(t(beta)%*%t(Zmat)))
      Y <- as.numeric(CenTime>=t)
      S0 <- as.numeric(t(Y)%*%M)
      YM <- Y*M
      S1 <- as.vector(t(YM)%*%Zmat)
      S2 <- t(Zmat)%*%diag(YM)%*%Zmat
      I <- I + K*(S2/S0-S1%*%t(S1)/(S0^2))
    }
    betaV[,i]<-as.vector(I/n)
    #print(eigen(I/n)$values)
  }
  return(list(betaMatrix=betaMatrix, betaV=betaV, countV=countV))
}

##
## @fn estBetaTime(h, data, timePoint, Nmin, coefTrue)
## @brief calculate the local constant partial maximum likelihood estimate with bandwidth 
##        h and fixed time point t.
## @return estimate of the coefficients at time point t
##
estBetaTime<-function(h, t, data, index, beta0.ini)
{
  beta0.ini.index <- beta0.ini[index]
  out.nlm <- tryCatch( nlm(loclikelihood1, beta0.ini.index,
                          t=t, h=h, Zmat=data$Zmat, CenTime=data$CenTime, Delta=data$Delta), 
                          error=function(e) list(estimate=rep(NA,p)))
  beta <- out.nlm$estimate
  return(beta)
}

##
## @fn loclikelihood1(beta, t, h, Zmat, CenTime, Delta)
## @brief calculate partial log likelihood of Cox model with estimated regression
##        coefficients 
## @return numeric, partial log likelihood of Cox model
##
loglikeCoxEst = function(Coef, timePoint, data_t, lbound, rbound)
{
  CenTime <- data_t$CenTime
  Delta <- data_t$Delta
  Zmat <- data_t$Zmat
  loglike <- 0
  for (i in which(CenTime<=rbound & CenTime>=lbound & Delta==1))
  {
    t <- CenTime[i]
    t1 <- max(timePoint[timePoint<=t]) 
    t2 <- min(timePoint[timePoint>=t])
    beta <- as.vector((Coef[,timePoint==t1]+Coef[,timePoint==t2])/2)
    M <- exp(as.vector(t(beta)%*%t(Zmat))) 
    Y <- as.numeric(CenTime>=t)
    G <- as.numeric(t(Y)%*%M)    
    loglike <- loglike+as.numeric(t(beta)%*%Zmat[i,]-log(G))  
  }
  return(loglike)
}

##
## @fn PE(h, kappa, data, fun, index)
## @brief calculate the prediction error (PE) by kappa-fold cross-validation method with 
##        a fix bandwidth h [Tian and Zucker:03]. 
## @return numeric, average prediction error of Cox model
##
PE = function(h, kappa, data, fun, index, beta0.ini)
{
  # keep only covariates in index set
  data.index <- data
  data.index$Zmat <- data$Zmat[,index]
  n <- nrow(data.index$Zmat)
  k <- floor(n/kappa)
  INDEX <- lapply(1:kappa, function(i, k){((i-1)*k+1):(i*k)}, k=k)
  INDEX[[kappa]] <- ((kappa-1)*k+1):n
  PEvec <- rep(NA, kappa)
  for (i in 1:kappa)
  {
    #data_est is the data use to estimate beta
    data_est <- data.index
    data_est$Zmat <- data.index$Zmat[-INDEX[[i]],]
    data_est$Delta <- data.index$Delta[-INDEX[[i]]]
    data_est$CenTime <- data.index$CenTime[-INDEX[[i]]]
    #data_test is the data use for test
    data_test <- data.index
    data_test$Zmat <- data.index$Zmat[INDEX[[i]],]
    data_test$Delta <- data.index$Delta[INDEX[[i]]]
    data_test$CenTime <- data.index$CenTime[INDEX[[i]]]
    #Calculate prediction error PE
    PE <- 0
    for (j in which(data_test$Delta==1))
    {
      tp <- data_test$CenTime[j]
      beta_s <- fun(h, tp, data_est, index, beta0.ini)
      M <- exp(as.vector(t(beta_s)%*%t(data_test$Zmat)))
      Y <- as.numeric(data_test$CenTime>=tp)
      PE <- PE-as.numeric(t(beta_s)%*%data_test$Zmat[j,]-log(as.numeric(t(Y)%*%M)))
    }
    PEvec[i] <- PE
  }
  return(mean(PEvec))
}

##
## @fn PE(h, kappa, data, fun, index)
## @brief search optimal bandwidth from a grid values of h by minimizing prediction error. 
## @return list with the following key-->value pairs:
##         bandwidthVect-->grid values of bandwidth h
##         APE-->vector of average prediction errors
##
searchBandwidth = function(kappa, data, bandwidthVect, fun, index, beta0.ini)
{
  APE <- lapply(bandwidthVect, PE, kappa=kappa, data=data, fun=fun, index=index, 
                beta0.ini=beta0.ini)
  return(list(bandwidthVect=bandwidthVect, APE=APE))
}

##
## @fn iniCI(ini.out, data, coefMatrix, timePoint, rbound, lbound, n, p, h1, Mpermu)
## @brief calculate the variance covariance matrix of initial estimator 
## @return list with the following key-->value pairs:
##
iniCI<-function(ini.out, data, coefMatrix, timePoint, rbound, lbound, n, p, h1, Mpermu){
  ntimes<-length(timePoint)
  # calculate the variance covariance matrix and corresponding 95% CI
  betaV<-ini.out$betaV
  # inverse Sigma(t)
  invbetaV<-matrix(NA, nrow(betaV), ncol(betaV))
  # variance-covariance matrix of m tilde
  sigmaM<-matrix(0, nrow=p, ncol=p)
  # standard error of beta tilde
  SEBeta<-matrix(NA, nrow=p, ncol=ntimes)
  #print(ntimes)
  for (k in 1:ntimes)
  {
    # use Moore-Penrose inverse to take care possible singular issue
    inv_sigma<-Ginv(matrix(betaV[,k], nrow=p, ncol=p))
    # inv_sigma<-solve(matrix(betaV[,k], nrow=p, ncol=p))
    invbetaV[,k]<-as.vector(inv_sigma)
    sigmaM<-sigmaM+inv_sigma
    # int K^2 from -1 to 1 =3/5
    #print(diag(inv_sigma))
    diag_inv_sigma<-pmax(diag(inv_sigma), 0)
    SEBeta[,k]<-sqrt(diag_inv_sigma*3/5/(n*h1))
    # print(diag(inv_sigma))
  }
  sigmaM<-sigmaM/(rbound-lbound)/ntimes/n
  # standard error of mtilde
  SEM<-sqrt(diag(sigmaM))
  
  # pointwise CI of betatilde
  ini.KGNG.ub.ci<-coefMatrix+1.96*SEBeta
  ini.KGNG.lb.ci<-coefMatrix-1.96*SEBeta
  # Confidence band
  tildeUM<-as.list(1:ntimes)
  for(i in 1:ntimes){
    numUtilde<-length(which(data$Delta==1 & abs((data$CenTime-timePoint[i])/h1)<1))
    tildeUM[[i]]<-matrix(NA, p, numUtilde)
    beta_t<-coefMatrix[,i]
    iter<-0
    for(j in which(data$Delta==1 & abs((data$CenTime-timePoint[i])/h1)<1)){
      iter<-iter+1
      K<-EpanKernel((data$CenTime[j]-timePoint[i])/h1)/h1
      M<-exp(as.vector(t(beta_t)%*%t(data$Zmat)))
      Y<-as.numeric(data$CenTime>=data$CenTime[j])
      S0<-as.numeric(t(Y)%*%M)
      YM<-Y*M
      S1<-as.vector(t(YM)%*%data$Zmat)
      tildeUM[[i]][,iter]<-(n)^(-1/2)*(h1)^(1/2)*(data$Zmat[j,]-S1/S0)*K
    }
  }
  # SEBeta is the inverse hat w(t)
  Gmatrix<-matrix(rnorm(Mpermu*n, 0, 1), Mpermu, n)
  Stilde<-matrix(NA, p, Mpermu)
  for (k in 1:Mpermu){
    G<-Gmatrix[k,]
    StildeM<-matrix(NA, p, ntimes)
    for(i in 1:ntimes){
      w_t<-1/SEBeta[,i]
      tildeU<-as.vector(tildeUM[[i]]%*%
                          G[which(data$Delta==1 & abs((data$CenTime-timePoint[i])/h1)<1)])
      StildeM[,i]<-abs(w_t*as.vector(matrix(invbetaV[,i], nrow=p, ncol=p)%*%tildeU))/
        (n*h1)^(1/2)
    }
    Stilde[,k]<-apply(StildeM, 1, max)
  }
  c_alpha<-apply(Stilde, 1, quantile, probs=0.90)
  ini.KGNG.ub.cb<-coefMatrix+diag(c_alpha)%*%SEBeta
  ini.KGNG.lb.cb<-coefMatrix-diag(c_alpha)%*%SEBeta
  
  return(list(ini.KGNG.ub.ci=ini.KGNG.ub.ci, ini.KGNG.lb.ci=ini.KGNG.lb.ci, 
              ini.KGNG.ub.cb=ini.KGNG.ub.cb, ini.KGNG.lb.cb=ini.KGNG.lb.cb,
              SEBeta=SEBeta, SEM=SEM, c_alpha=c_alpha, invbetaV=invbetaV, sigmaM=sigmaM)) 
}

##
## @fn prelimStep(data, n, p, beta0.ini, nlength, timePoint, Nmin, h, lbound, rbound, thetaVect)
## @brief implement preliminary step 
##
prelimStep <- function(data, n, p, beta0.ini, nlength, timePoint, Nmin, h, 
                        lbound, rbound, thetaVect, positive){   
  ptm <- proc.time()
  # estimate coefficients at first step
  coefini =  beta0.ini%*%t(rep(1, nlength))
  coefMatrix <- estBeta_h(h, data, timePoint, Nmin, coefini)$betaMatrix   
  # selection between zero and non-zero covariates 
  # get Y and X matrix of Cox model with time-dependent covariates
  k<-0
  for(i in 1:n)
  {
    for (j in 1:(nlength-1))
    {
      if (timePoint[j]<data$CenTime[i]) k <- k+1  
    }
  }  
  start <- rep(NA, k)
  stop <- rep(NA, k)
  cen.time <- rep(NA, k)
  X <- matrix(NA, k, p)  
  k2 <- 0
  for(i in 1:n)
  {
    for (j in 1:(nlength-1))
    {
      if (timePoint[j]<data$CenTime[i]) 
      {
        k2 <- k2+1
        start[k2] <- timePoint[j]
        stop[k2] <- timePoint[(j+1)]
        if(timePoint[(j+1)]>=data$CenTime[i]) cen.time[k2] = data$Delta[i] else 
          cen.time[k2] = 0
        X[k2,1:p] <- data$Zmat[i,]*(coefMatrix[,j]+coefMatrix[,(j+1)])/2
      }  
    }
  }  
  Y <- Surv(start,stop,cen.time) #Define Y as a survival subject 
  # tuning with BIC to find optimal theta*
  bic.OPT <- NA
  startbeta <- rep(0,p)
  for(i in 1:length(thetaVect))
  {
    # use "penalized" to get lambda_hat
    pen <- penalized(Y, penalized = X, lambda1=thetaVect[i], positive=positive, 
                     startbeta=startbeta) 
    lambda_fixed <- coefficients(pen,"penalized")
    startbeta = lambda_fixed
    Index.nonzero.temp <- which(lambda_fixed!=0)                  
    # calculate BIC value
    logLike <- loglik(pen)
    df <- sum(lambda_fixed!=0)
    bic <- -2*logLike+df*log(n*h/(rbound-lbound))/(h/(rbound-lbound))      
    # estimated coefficients 
    coefMatrix2_temp <- diag(lambda_fixed)%*%coefMatrix  
    if (i==1){
      bic.OPT <- bic
      Index.nonzero <- Index.nonzero.temp
      coefMatrix2 <- coefMatrix2_temp
      theta.star.opt = thetaVect[i]
    } else if (bic<bic.OPT){
      bic.OPT<-bic
      Index.nonzero <- Index.nonzero.temp
      coefMatrix2 <- coefMatrix2_temp 
      theta.star.opt = thetaVect[i]
    }
  }
  # Stop the clock
  time.prelim = print(proc.time() - ptm)
  ptm <- proc.time() 
  return(list(Index.nonzero=Index.nonzero, coefMatrix=coefMatrix, coefMatrix2=coefMatrix2,
              theta.star.opt=theta.star.opt, time.prelim=time.prelim))
}