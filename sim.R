####### Project 3 - Direct maximization
####### Final R Code - Simulation Functions ########
####### Date: 3/25/19

#### Purpose of Code: Simulations using method. First set of simulations used bivarate normal
###                   data

library(Rsolnp)
library(mvtnorm)
library(aucm)

###### A:  Simulation Functions ########
simOne <- function(n,CB,intVal=NULL,hPower=0.5,prev=NULL,ftol=1e-4,delta=1e-8,simOne.seed,
                   datFun,prevCalc=TRUE){
  trainDat <- datFun(n=n, prev = prev, seed=simOne.seed )
  trainY <- trainDat$y
  trainX <- trainDat$x
  trainXD <- trainDat$x[trainDat$y==1,]
  trainXND <- trainDat$x[trainDat$y==0,]
  
  testSeed <- simOne.seed + 112312
  testDat <- datFun(n=100000, prev = prev, seed=testSeed )
  testY <- testDat$y
  testX <- testDat$x
  testXD <- testDat$x[testDat$y==1,]
  testXND <- testDat$x[testDat$y==0,]
  
  if(prevCalc==TRUE){
    prev <- mean(trainY)
  }
  
  w <- CB * ((1-prev)/prev)
  ### Implemented Comparater Methods ###
  ## Logistic Regression
  log <- glm(trainY ~ trainX, family = "binomial")
  log.coef <- log$coef
  # Optimized t 
  stdLogit <- logistOpt(coef = log.coef,w = w,
                        y = trainY, x = trainX,CB = CB,hPower = hPower,prev = prev)
  std.pHat.D <- expit(cbind(1,testXD) %*% log.coef)
  std.pHat.ND <- expit(cbind(1,testXND) %*% log.coef)
  
  sNB.log.test <- ind.sNB.thresh(par = stdLogit["t.opt"],
                                 xD = std.pHat.D,xND = std.pHat.ND,w = w)
  tpr.log.test <- mean(std.pHat.D > stdLogit["t.opt"])
  fpr.log.test <- mean(std.pHat.ND > stdLogit["t.opt"])
  # Using cost benefit ratio
  sNB.logCB.test <- ind.sNB.thresh(par = stdLogit["t.cb"],
                                   xD = std.pHat.D,xND = std.pHat.ND,w = w)
  tpr.logCB.test <- mean(std.pHat.D > stdLogit["t.cb"])
  fpr.logCB.test <- mean(std.pHat.ND > stdLogit["t.cb"])
  
  stdRes <- list("stdOpt" = stdLogit, 
                 "sNB.log.test" = sNB.log.test,
                 "tpr.log.test" = tpr.log.test, "fpr.log.test" = fpr.log.test,
                 "sNB.logCB.test" = sNB.logCB.test,
                 "tpr.logCB.test" = tpr.logCB.test, "fpr.logCB.test" = fpr.logCB.test)
  
  
  ## robust logistic regression estimates
  robustLog <- rlogit(trainY ~ trainX,dat = trainDat)
  rLog.coef <- robustLog$coef
  # Optimized t 
  robLogit <- logistOpt(coef = rLog.coef, w = w,
                        y = trainY, x = trainX,CB = CB,hPower = hPower,prev = prev)
  rob.pHat.D <- expit(cbind(1,testXD) %*% rLog.coef)
  rob.pHat.ND <- expit(cbind(1,testXND) %*% rLog.coef)
  
  sNB.roblog.test <- ind.sNB.thresh(par = robLogit["t.opt"],
                                    xD = rob.pHat.D,xND = rob.pHat.ND,w = w)
  tpr.roblog.test <- mean(rob.pHat.D > robLogit["t.opt"])
  fpr.roblog.test <- mean(rob.pHat.ND > robLogit["t.opt"])
  # Using cost benefit ratio
  sNB.roblogCB.test <- ind.sNB.thresh(par = robLogit["t.cb"],
                                      xD = rob.pHat.D,xND = rob.pHat.ND,w = w)
  tpr.roblogCB.test <- mean(rob.pHat.D > robLogit["t.cb"])
  fpr.roblogCB.test <- mean(rob.pHat.ND > robLogit["t.cb"])
  
  robRes <- list("robOpt" = robLogit, 
                 "sNB.log.test" = sNB.roblog.test,
                 "tpr.log.test" = tpr.roblog.test, "fpr.log.test" = fpr.roblog.test,
                 "sNB.logCB.test" = sNB.roblogCB.test,
                 "tpr.logCB.test" = tpr.roblogCB.test, "fpr.logCB.test" = fpr.roblogCB.test)
  
  ### Direct optimization 
  doMax <- maxSNB(dat = trainDat,y = trainY,x = trainX, 
                  CB = CB,intVal = intVal,hPower = hPower,prev = prev,
                  ftol = ftol,delta = delta)
  
  sNB.test <- ind.SNB.opt(theta = doMax$directOpt.max$theta,xD = testXD,xND = testXND,w = w) 
  tpr.test <- tpr.opt(theta= doMax$directOpt.max$theta,xD = testXD,xND = testXND)
  fpr.test <- fpr.opt(theta= doMax$directOpt.max$theta,xD = testXD,xND = testXND)
  doRes <- list("doMax" = doMax, 
                "sNB.test" = sNB.test,
                "tpr.test" = tpr.test, "fpr.test" = fpr.test)
  
  
  #### Direct optimization with updated optimization 
  UsNB.test <- ind.SNB.opt(theta = doMax$directOpt.max$uRes$uTheta,xD = testXD,xND = testXND,w = w) 
  Utpr.test <- tpr.opt(theta= doMax$directOpt.max$uRes$uTheta,xD = testXD,xND = testXND)
  Ufpr.test <- fpr.opt(theta= doMax$directOpt.max$uRes$uTheta,xD = testXD,xND = testXND)
  UdoRes <- list("UdoMax" = doMax$directOpt.max$uRes, 
                 "sNB.test" = UsNB.test,
                 "tpr.test" = Utpr.test, "fpr.test" = Ufpr.test)
  
  
  res <- list("stdRes" = stdRes, "robRes" = robRes, "doRes" = doRes, "UdoRes" = UdoRes)
  return(res)
}

simMany <- function(B,n,CB,intVal=NULL,hPower,prev=NULL,prevCalc=TRUE,
                    ftol=1e-8,delta=1e-7,int.seed,
                    datFun,thetaLength){
  stdLR <- data.frame(matrix(NA,nrow=B,ncol=15))
  names(stdLR) <- c("t.opt","t.opt.conv","t.cb", "sNB.log.train","tpr.log.train",
                    "fpr.log.train","sNB.logCB.train","tpr.logCB.train","fpr.logCB.train",
                    "sNB.log.test","tpr.log.test",
                    "fpr.log.test","sNB.logCB.test","tpr.logCB.test","fpr.logCB.test")
  
  robLR <- data.frame(matrix(NA,nrow=B,ncol=15))
  names(robLR) <- c("t.opt","t.opt.conv","t.cb", "sNB.log.train","tpr.log.train",
                    "fpr.log.train","sNB.logCB.train","tpr.logCB.train","fpr.logCB.train",
                    "sNB.log.test","tpr.log.test",
                    "fpr.log.test","sNB.logCB.test","tpr.logCB.test","fpr.logCB.test")
  
  nColInt <- thetaLength + 2
  doOptInt <- data.frame(matrix(as.numeric(NA),nrow=B,ncol=nColInt))
  names(doOptInt)[(thetaLength+1):nColInt] <- c("t0.conv","stdInt") #names for all but theta entries
  
  nColOpt <- thetaLength + 8
  doOpt <- data.frame(matrix(as.numeric(NA),nrow=B,ncol=nColOpt))
  names(doOpt)[(thetaLength+1):nColOpt] <- c("obj","conv","sNB.train","tpr.train","fpr.train","sNB.test","tpr.test","fpr.test")
  
  UnColOpt <- 8
  UdoOpt <- data.frame(matrix(as.numeric(NA),nrow=B,ncol=UnColOpt))
  names(UdoOpt) <- c("uT","uT.conv","UsNB.train","Utpr.train","Ufpr.train",
                     "UsNB.test","Utpr.test","Ufpr.test")
  
  for(i in 1:B){
    simMany.seed <- int.seed + i 
    oneRep <- simOne(n = n, CB = CB,intVal = intVal,
                     hPower = hPower,prev = prev,
                     ftol = ftol,delta = delta,
                     simOne.seed = simMany.seed,datFun = datFun,prevCalc = prevCalc)
    stdLR[i,] <- c(oneRep$stdRes$stdOpt,
                   oneRep$stdRes$sNB.log.test,oneRep$stdRes$tpr.log.test,oneRep$stdRes$fpr.log.test,
                   oneRep$stdRes$sNB.logCB.test,oneRep$stdRes$tpr.logCB.test,oneRep$stdRes$fpr.logCB.test)
    
    robLR[i,] <- c(oneRep$robRes$robOpt,
                   oneRep$robRes$sNB.log.test,oneRep$robRes$tpr.log.test,oneRep$robRes$fpr.log.test,
                   oneRep$robRes$sNB.logCB.test,oneRep$robRes$tpr.logCB.test,oneRep$robRes$fpr.logCB.test)
    
    doOptInt[i,] <- c(oneRep$doRes$doMax$intVal,oneRep$doRes$doMax$t0.conv,NA)
    doOptInt$stdInt <- ifelse(oneRep$doRes$doMax$whichMax=="std",1,0)
    
    doOpt[i,] <- c(oneRep$doRes$doMax$directOpt.max$theta,
                   oneRep$doRes$doMax$directOpt.max$obj,oneRep$doRes$doMax$directOpt.max$conv,
                   oneRep$doRes$doMax$directOpt.max$trainRes,
                   oneRep$doRes$sNB.test,oneRep$doRes$tpr.test,oneRep$doRes$fpr.test)
    
    UdoOpt[i,] <- c(oneRep$UdoRes$UdoMax$uT,oneRep$UdoRes$UdoMax$uTheta.conv,
                    oneRep$UdoRes$UdoMax$UtrainRes,
                    oneRep$UdoRes$sNB.test,oneRep$UdoRes$tpr.test,oneRep$UdoRes$fpr.test)
    print(i)
  }
  names(doOptInt)[1:thetaLength] <- names(oneRep$doRes$doMax$intVal)
  names(doOpt)[1:thetaLength] <- names(oneRep$doRes$doMax$directOpt.max$theta)
  
  ## Summarizes results 
  #comparator results
  stdLR.mean <- apply(stdLR[,-2],2,mean,na.rm = T)
  stdLR.sd <- apply(stdLR[,-2],2,sd,na.rm = T)
  stdLR.prop <- mean(stdLR[,2]==0)
  
  robLR.mean <- apply(robLR[,-2],2,mean,na.rm = T)
  robLR.sd <- apply(robLR[,-2],2,sd,na.rm = T)
  robLR.prop <- mean(robLR[,2]==0)
  
  compMethod.mean <- data.frame(rbind(stdLR.mean,robLR.mean))
  compMethod.sd <- data.frame(rbind(stdLR.sd,robLR.sd))
  names(compMethod.mean) <- names(compMethod.sd) <- names(robLR)[-2]
  compMethod.t0prop <- c(stdLR.prop,robLR.prop)
  names(compMethod.t0prop) <- c("std.t0.propConv", "rob.t0.propConv")
  compRes <- list("compMethod.mean" = compMethod.mean,
                  "compMethod.sd" = compMethod.sd,
                  "compMethod.t0prop" = compMethod.t0prop)
  
  #direct opt results
  doOptInt.mean <- apply(doOptInt[,1:thetaLength],2,FUN = mean,na.rm = T)
  doOptInt.sd <- apply(doOptInt[,1:thetaLength],2,FUN = sd,na.rm = T)
  doOptInt.t0prop <- mean(doOptInt$t0.conv==0)
  doOptInt.stdInt <- mean(doOptInt$stdInt)
  
  doOpt.mean <- apply(doOpt[,!(names(doOpt) %in% "conv")],2,mean,na.rm = T)
  doOpt.sd <- apply(doOpt[,!(names(doOpt) %in% "conv")],2,sd,na.rm = T)
  doOpt.conv <- mean(doOpt$conv==0)
  
  #direct updated  results
  UdoOpt.mean <- apply(UdoOpt[,!(names(UdoOpt) %in% "uT.conv")],2,mean,na.rm = T)
  UdoOpt.sd <- apply(UdoOpt[,!(names(UdoOpt) %in% "uT.conv")],2,sd,na.rm = T)
  UdoOpt.conv <- mean(UdoOpt$uT.conv==1)
  
  doRes <- list("doOptInt.mean" = doOptInt.mean, "doOptInt.sd" = doOptInt.sd, 
                "doOptInt.t0prop" = doOptInt.t0prop, "doOptInt.stdInt" = doOptInt.stdInt,
                "doOpt.mean" = doOpt.mean, "doOpt.sd" = doOpt.sd, "doOpt.conv" = doOpt.conv,
                "UdoOpt.mean" = UdoOpt.mean, "UdoOpt.sd" = UdoOpt.sd, "UdoOpt.conv" = UdoOpt.conv)  
  
  res <- list("doRes" = doRes, "compRes" = compRes,
              "doOptFull" = doOpt, "doOptInt" = doOptInt, "stdFull" = stdLR, "robFull" = robLR,
              "UdoOptFull" = UdoOpt)
  return(res)
}



###### 1: Simulation for Std Bivariate Normal Case #######
stdBiv <- biVar.simData <- function(n,prev,mu.D = c(1,1),sigmaD = diag(1,nrow = 2),
                                    sigmaND = diag(1,nrow=2),seed){
  set.seed(seed)
  
  #first simulate outcome 
  d <- rbinom(n = n,size = 1,prob = prev)
  d <- sort(d)
  #now simulate bivariate markers
  x.D <- rmvnorm(n = sum(d==1),mean = mu.D,sigma = sigmaD)
  x.ND <- rmvnorm(n = sum(d==0),mean = c(0,0),sigma = sigmaND)
  
  dat <- data.frame("y" = d, "x1"= c(x.ND[,1],x.D[,1]), "x2"= c(x.ND[,2],x.D[,2])) 
  res <- list("dat" = dat, "y" = d, "x" = cbind(c(x.ND[,1],x.D[,1]),c(x.ND[,2],x.D[,2])) )
  return(res)
}


#### N = 500
### This line will run a single set of simulations 
stdNorm.n500.w1 <- simMany(B = 500, n = 500,CB = (1/9),
                           intVal = NULL,hPower = 0.5,prev = 0.1,int.seed = 349,
                           datFun = stdBiv,thetaLength = 3)

