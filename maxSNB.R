####### Project 3 - Direct maximization
####### Final R Code - Functions ########
####### Date: 3/25/19


#### Purpose of Code: Functions used for direct optimization method 
###                   Run this whole file before running simulations

library(Rsolnp)
library(mvtnorm)
library(aucm)

########## A: Functions Used for Optimizations ############
logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}

### Optimization Functions for sNB
## smooth approximation function of sNB 
smooth.SNB.opt <- function(theta,xD,xND,sn,w){
  zD <- as.matrix(xD) %*% as.matrix(theta[-length(theta)])
  zND <- as.matrix(xND) %*% as.matrix(theta[-length(theta)])
  
  sTPR <- mean(pnorm(q = (zD-theta[length(theta)])/sn))
  sFPR <- mean(pnorm(q = (zND-theta[length(theta)])/sn))
  sSNB <- -1*(sTPR -  w*sFPR)  
  return(sSNB)
}


# constraint function used in optimzation
const <- function(theta,xD,xND,sn,w){
  res <- sqrt(sum(theta[-length(theta)]^2))
  return(res)
}

# t
tpr.opt <- function(theta,xD,xND){
  z.d <- as.matrix(xD) %*% as.matrix(theta[-length(theta)])
  z.nd <-   as.matrix(xND) %*% as.matrix(theta[-length(theta)])
  
  TPR <- mean(z.d > theta[length(theta)])
  return(TPR)
}

# empirical estimator of sNB, TPR and FPR
ind.SNB.opt <- function(theta,xD,xND,sn,w){
  z.d <- as.matrix(xD) %*% as.matrix(theta[-length(theta)])
  z.nd <-   as.matrix(xND) %*% as.matrix(theta[-length(theta)])
  
  TPR <- mean(z.d > theta[length(theta)])
  FPR <- mean(z.nd > theta[length(theta)])
  SNB <- (TPR -  w*FPR)  
  return(SNB)
}

fpr.opt <- function(theta,xD,xND,sn,w){
  z.d <- as.matrix(xD) %*% as.matrix(theta[-length(theta)])
  z.nd <-   as.matrix(xND) %*% as.matrix(theta[-length(theta)])
  
  FPR <- mean(z.nd > theta[length(theta)])
  return(FPR)
}

# plot of NB for varying risk threshold
sNBplot <- function(tVec,xD,xND,w,sn){
  sNBvec <- cbind(tVec,NA)
  for(i in 1:length(tVec)){
    sNBvec[i,2] <- smooth.SNB.thresh(par = tVec[i],xD = xD,xND = xND,w = w,sn = sn)
  }
  return(sNBvec)
}

## Optimization functions for risk threshold 
smooth.SNB.thresh <- function(par,xD,xND,w,sn){
  sTPR <- mean(pnorm(q = (xD-par)/sn))
  sFPR <- mean(pnorm(q = (xND-par)/sn))
  sSNB <- (sTPR -  (w*sFPR))
  return(sSNB)
}
ind.sNB.thresh <- function(par,xD,xND,w){
  TPR <- mean(xD > par)
  FPR <- mean(xND > par)
  sNB <- (TPR -  (w*FPR))
  return(sNB)
}

##### B: Optimization Routines - Direct Opt #######
directOpt <- function(y,x,int,w,CB,hPower,prev,ftol,delta){
  xD <- x[y==1,]
  xND <- x[y==0,]
  
  sn <- sd(x %*% int[-length(int)])/nrow(x)^hPower
  optim.smooth <- solnp(pars = int,fun = smooth.SNB.opt,eqfun = const,eqB = 1,
                        xD = xD, xND = xND, sn = sn, w=w,
                        control = list(trace=0,tol=ftol,delta=delta))
  theta <- optim.smooth$pars
  names(theta)[length(theta)] <- "t"
  ## sNB under my method 
  sNB.train <- ind.SNB.opt(theta = optim.smooth$pars,xD = xD,xND = xND,w = w) 
  tpr.train <- tpr.opt(theta= optim.smooth$pars,xD = xD,xND = xND)
  fpr.train <- fpr.opt(theta= optim.smooth$pars,xD = xD,xND = xND)
  trainRes <- c(sNB.train,tpr.train,fpr.train)
  names(trainRes) <- c("sNB.train","tpr.train","fpr.train")
  
  ### What if I did another optimization to find t after giving theta estimates
  snUpdate <- sd(x %*% theta[-length(theta)])/nrow(x)^hPower
  xD.do <- xD %*% theta[-length(theta)]
  xND.do <- xND %*% theta[-length(theta)]
  updateT <- optimize(f = smooth.SNB.thresh,
                      interval = c(min(c(xD.do,xND.do)),max(c(xD.do,xND.do))),
                      xD = xD.do,xND = xND.do,w = w,sn = snUpdate,maximum = TRUE)
  uTheta <- c(theta[-length(theta)],updateT$maximum)
  uT <- uTheta[length(uTheta)]
  sNB.train <- ind.SNB.opt(theta = uTheta,xD = xD,xND = xND,w = w) 
  tpr.train <- tpr.opt(theta= uTheta,xD = xD,xND = xND)
  fpr.train <- fpr.opt(theta= uTheta,xD = xD,xND = xND)
  UtrainRes <- c(sNB.train,tpr.train,fpr.train)
  names(UtrainRes) <- c("UsNB.train","Utpr.train","Ufpr.train")
  uRes <- list("uTheta" = uTheta, "uT" = uT, 
               "uTheta.conv" = ifelse(uT < 0,0,1),
               "UtrainRes" = UtrainRes)
  
  res <- list("theta"= theta,"obj" = optim.smooth$values[length(optim.smooth$values)],
              "conv" = optim.smooth$convergence,
              "trainRes" = trainRes,"uRes" = uRes)
  return(res)
}

logistOpt <- function(coef,y,x,w,CB,hPower,prev){
  xD <- x[y==1,]
  xND <- x[y==0,]
  
  #converting to risk scale
  zD <- as.matrix(cbind(1,xD)) %*% coef 
  pHat.D  <- expit(zD)
  zND <-  as.matrix(cbind(1,xND)) %*% coef 
  pHat.ND <- expit(zND)
  sn <- sd(c(zD,zND))/((length(c(zD,zND)))^(hPower))
  
  ## getting omega
  log.t.cb <- CB/(1+CB)
  
  
  ## finding optimized t given logistic parameters
  log.t.opt <- optim(par = log.t.cb,fn = smooth.SNB.thresh,method = "L-BFGS-B",
                     control=list(fnscale=-1),
                     xD = pHat.D, xND = pHat.ND,
                     w = w,lower=0,upper = 1,sn=sn)
  t.opt <- log.t.opt$par
  if(log.t.opt$convergence!=0){
    log.t.optGrid <- cbind(seq(0,1,0.001),
                           sapply(X = seq(0,1,0.001),smooth.SNB.thresh,
                                  xD = pHat.D,xND = pHat.ND,w = w,sn = sn))
    t.opt <- log.t.optGrid[which.max(log.t.optGrid[,2]),1]
  }
  
  sNB.log.train <- ind.sNB.thresh(par = t.opt,xD = pHat.D,xND = pHat.ND,w = w)
  tpr.log.train <- mean(pHat.D > t.opt)
  fpr.log.train <- mean(pHat.ND > t.opt)
  
  
  sNB.logCB.train <- ind.sNB.thresh(par = log.t.cb,xD = pHat.D,xND = pHat.ND,w = w)
  tpr.logCB.train <- mean(pHat.D > log.t.cb)
  fpr.logCB.train <- mean(pHat.ND > log.t.cb)
  
  res <- c(t.opt,log.t.opt$convergence, log.t.cb,
           sNB.log.train, tpr.log.train, fpr.log.train,
           sNB.logCB.train, tpr.logCB.train, fpr.logCB.train)
  names(res) <- c("t.opt","t.opt.conv", "t.cb",
                  "sNB.log.train","tpr.log.train","fpr.log.train", 
                  "sNB.logCB.train","tpr.logCB.train", "fpr.logCB.train")
  return(res)
}

maxSNB <- function(dat,y,x,CB,intVal=NULL,hPower,prev,ftol=1e-6,delta=1e-8){
  xD <- x[y==1,]
  xND <- x[y==0,]
  
  w <- CB * ((1-prev)/prev)
  #### (1) Getting inital values 
  ## logistic coefficient estimates under logistic regression
  log <- glm(y ~ x, family = "binomial")
  log.coef <- log$coef
  stdOpt <- logistOpt(coef = log.coef,y = y,x = x,CB = CB,hPower = 0.5,prev = prev, w = w)
  
  ## robust logistic regression estimates
  robustLog <- rlogit(y ~ x,dat = dat)
  rLog.conv <- ifelse(robustLog$convergence=="TRUE",1,0)
  rLog.coef <- robustLog$coef
  robOpt <- logistOpt(coef = rLog.coef,y = y,x = x,CB = CB,hPower = 0.5,prev = prev, w = w)
  
  
  if(is.null(intVal)){
    if(robOpt["sNB.log.train"] > stdOpt["sNB.log.train"]){
      rLog.coefStd  <- rLog.coef[-1]/norm(rLog.coef[-1],type="2")
      xD.robLog <- xD %*% rLog.coefStd
      xND.robLog <- xND %*% rLog.coefStd
      snMax <- sd(c(xND.robLog,xD.robLog))/((nrow(x))^(hPower))
      
      rLogStd.t.opt <- optimize(f = smooth.SNB.thresh,
                                xD = xD.robLog, xND = xND.robLog, w = w,sn = snMax,
                                interval = c(min(c(xD.robLog,xND.robLog)),max(c(xD.robLog,xND.robLog))),
                                maximum = TRUE)
      
      t0 <- rLogStd.t.opt$maximum
      t0.conv <- ifelse(rLogStd.t.opt$objective <= 0, 0, 1)
      #if(t0.conv!=1){
      #  tSeq <- seq(min(c(xD.robLog,xND.robLog)),max(c(xD.robLog,xND.robLog)),0.001)
      #  t0.grid <- cbind(tSeq,sapply(X = tSeq,smooth.SNB.thresh,
      #                               xD = xD.robLog,xND = xND.robLog,w = w,sn = snMax))
      #  t0 <- t0.grid[which.max(t0.grid[,2]),1]
      #}
      maxCoef <-  c(rLog.coefStd,t0)
      whichCoef <- "rob"
    }else{ 
      log.coefStd  <- log.coef[-1]/norm(log.coef[-1],type="2") 
      xD.stdLog <- xD %*% log.coefStd
      xND.stdLog <- xND %*% log.coefStd
      snMax <- sd(c(xD.stdLog,xND.stdLog))/((nrow(x))^(hPower))
      
      stdLog.t.opt <- optimize(f = smooth.SNB.thresh,
                               xD = xD.stdLog, xND = xND.stdLog, w = w,sn = snMax,
                               interval = c(min(c(xD.stdLog,xND.stdLog)),max(c(xD.stdLog,xND.stdLog))),
                               maximum = TRUE)
      
      t0 <- stdLog.t.opt$maximum
      t0.conv <- ifelse(stdLog.t.opt$objective <= 0, 0, 1)
      #if(stdLog.t.opt$convergence!=0){
      #  tSeq <- seq(min(c(xD.stdLog,xND.stdLog)),max(c(xD.stdLog,xND.stdLog)),0.001)
      #  t0.grid <- cbind(tSeq,sapply(X = tSeq,smooth.SNB.thresh,
      #                               xD = xD.stdLog,xND = xND.stdLog,w = w,sn = snMax))
      #  t0 <- t0.grid[which.max(t0.grid[,2]),1]
      #}
      
      maxCoef <- c(log.coefStd,t0)
      whichCoef <- "std"
    }
    names(maxCoef)[length(maxCoef)] <- "t0"
  }
  else{maxCoef <- intVal
  whichCoef <- NULL
  t0.conv <- NULL
  }
  
  #### Applying Direct Opt 
  directOpt.max <- directOpt(y = y,x = x,int = maxCoef,CB = CB,w = w,
                             hPower = hPower,prev = prev,ftol = ftol,delta = delta)
  res <- list("directOpt.max" = directOpt.max,
              "intVal" = maxCoef,
              "whichMax" = whichCoef,
              "t0.conv" = t0.conv)
  
  return(res)
}
