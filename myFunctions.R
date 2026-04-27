cal.DS <- function(myVar,nQtr,UniDate){
  
  DS <- list()
  DS$plus0 <- DS$plus1 <- DS$plus2 <- DS$plus3 <- DS$plus4 <- matrix(NA,nQtr,9)
  H <-  c('plus0','plus1','plus2','plus3','plus4')
  nH <- length(H) 
  nFcs <- matrix(NA,nQtr,nH)
  
  for (t in 1:nQtr){
    qtr <- UniDate[t]
    myVar.qtr <- subset(myVar, Date %in% qtr)
    
    DS$plus0[t,] <- cal.DS.idv(myVar.qtr$plus0)
    DS$plus1[t,] <- cal.DS.idv(myVar.qtr$plus1)
    DS$plus2[t,] <- cal.DS.idv(myVar.qtr$plus2)
    DS$plus3[t,] <- cal.DS.idv(myVar.qtr$plus3)
    DS$plus4[t,] <- cal.DS.idv(myVar.qtr$plus4)
    
    nFcs[t,] <- apply(myVar.qtr[H], 2, function(x) sum(!is.na(x))) 
  }
  
  DS <- DS[H]
  
  DS.summary <- lapply(DS,colMeans,na.rm = T) %>% unlist() %>% matrix(nrow=nH,ncol=9,byrow = T)
  
  row.names(DS.summary) <- H
  colnames(DS.summary) <- c('Mean','SD','Skewness','Kurtosis','10% Q','25% Q','50% Q','75% Q','90% Q')
  colnames(nFcs) <- H
  
  list(DS.summary=DS.summary, nFcs=nFcs)
  
}

cal.DS.idv <- function(myVec){
  DS.idv <- matrix(NA,1,9)
  
  DS.idv[1] <- mean(myVec, na.rm = T)
  DS.idv[2] <- sd(myVec, na.rm = T)
  DS.idv[3] <- skewness(myVec, na.rm = TRUE)
  DS.idv[4] <- kurtosis(myVec, na.rm = TRUE)+3
  
  DS.idv[5:9] <- quantile(myVec, c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = T)
  
  DS.idv
}

normalize.density <- function(myDens){
  y <- myDens$y / trapz(x = myDens$x, y = myDens$y)
  # y[y==0] <- 2.225074e-308 #.Machine$double.xmin
  # y <- RegulariseByAlpha(y, x = myDens$x, alpha = 0.01)
  y
}

mydens2lqd <-  function(dens, dSup, lqdSup = seq(0, 1, length.out = length(dSup)), t0 = dSup[1], verbose = TRUE){
  
  # Check density requirements
  if(any(dens < 0)){
    stop('Please correct negative density values.')
  }
  
  if(abs( trapzRcpp(X = dSup, Y = dens) - 1) > 1e-5){
    
    warning('Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.')
    dens = dens/trapzRcpp(X = dSup, Y = dens)
    
  }
  
  if(any(dens == 0)){
    if(verbose){
      print("There are some zero density values - truncating support grid so all are positive")
    }
    lbd = min(which(dens > 0))
    ubd = max(which(dens > 0))
    dens = dens[lbd:ubd]
    dSup = dSup[lbd:ubd]
    dens = dens/trapzRcpp(X = dSup, Y = dens)
  }
  
  N = length(dSup)
  
  # Check LQD output grid
  if(is.null(lqdSup)){
    
    lqdSup = seq(0, 1, length.out = N)
    
  }else if(!all.equal( range(lqdSup),c(0,1) )){
    
    if(verbose){
      print("Problem with support of the LQD domain's boundaries - resetting to default.")
    }
    lqdSup = seq(0, 1, length.out = N)
    
  }
  
  # Check t0
  if(!(t0 %in% dSup)){
    
    if(verbose){
      print("t0 is not a value in dSup - resetting to closest value")
    }
    t0 = dSup[which.min(abs(dSup - t0))]
    
  }
  
  M = length(lqdSup) 
  c_ind = which(dSup == t0)
  
  # Get CDF and lqd on temporary grid, compute c
  tmp = cumtrapzRcpp(X = dSup, dens)
  c = tmp[c_ind]
  
  indL = duplicated(tmp[1:floor(N/2)])
  indR = duplicated(tmp[(floor(N/2)+1):N], fromLast = TRUE)
  qtemp = tmp[!c(indL, indR)]
  lqd_temp = -log(dens[!c(indL, indR)]);
  
  # Interpolate lqdSup, keeping track of Inf values at boundary, then compute c
  lqd = rep(0, 1, M)
  
  if(any(is.infinite(lqd_temp[c(1, N)]))){
    
    tmpInd = 1:N
    Ind = 1:M
    
    if(lqd_temp[1] == Inf){
      
      lqd[1] = Inf
      tmpInd = tmpInd[-1]
      Ind = Ind[-1]
      
    }
    
    if(lqd_temp[N] == Inf){
      
      lqd[M] = Inf
      tmpInd = tmpInd[-length(tmpInd)]
      Ind = Ind[-length(Ind)]
      
    }
    
    lqd[Ind] = approx(x = qtemp[tmpInd], y = lqd_temp[tmpInd], xout = lqdSup[Ind], rule = 2)$y 
    
  }else{
    
    lqd = approx(x = qtemp, y = lqd_temp, xout = lqdSup, rule = c(2,2))$y 
    
  }
  
  return(list('lqdSup',  lqdSup, 'lqd' = lqd, 'c' = c))
}

mylqd2dens <-  function(lqd, lqdSup = seq(0, 1, length.out = length(lqd)), dSup, t0 = 0, c = 0, useSplines = TRUE, cut = c(0, 0), verbose = TRUE){
  
  if(!all.equal( range(lqdSup),c(0,1) )){
    
    warning("Problem with support of the LQD domain's boundaries - resetting to default.")
    lqdSup = seq(0, 1, length.out = length(lqd))
    
  }
  
  M = length(lqd)
  r = which(exp(lqd) == Inf)
  
  if(length(r) > 0){
    
    if(any(r < floor(M/2))){
      cut[1] = max(cut[1], max(r[r < floor(M/2)]))
    }
    if(any(r >= floor(M/2))){
      cut[2] = max(cut[2], M - min(r[r >= floor(M/2)]) + 1)
    }    
  }
  
  # Cut boundaries
  lqdSup = lqdSup[(cut[1] + 1):(M - cut[2])]
  lqd = lqd[(cut[1] + 1):(M - cut[2])]
  M = length(lqd) # reset N
  
  if(!(c %in% lqdSup)){
    
    if(c < lqdSup[1] || c > lqdSup[M]){
      
      stop("c is not contained withing range of lqdSup after cutoff")
      
    }
    
    if(verbose){
      
      print("c is not equal to a value in lqdSup - resetting to closest value")
      
    }
    c = lqdSup[which.min(abs(lqdSup - c))]
    
  }
  
  c_ind = which(lqdSup == c)
  
  if( useSplines ){    # Could fit spline if this yields more accurate numerical integration
    
    lqd_sp = splinefun(lqdSup, lqd, method = 'natural')
    lqd_exp = function(t) exp(lqd_sp(t))
    # Get grid for density space
    dtemp = t0 + c(0, cumsum(sapply(2:length(lqdSup), function(i) integrate(lqd_exp, lqdSup[i - 1], lqdSup[i])$value))) - integrate(lqd_exp, lqdSup[1], lqdSup[c_ind])$value
    
  } else {
    # Get grid and function for density space
    dtemp = t0 + cumtrapzRcpp(lqdSup, exp(lqd)) - trapzRcpp(lqdSup[1:c_ind], exp(lqd[1:c_ind]))
  }
  
  # Remove duplicates
  indL = duplicated(dtemp[1:floor(M/2)], fromLast = TRUE)
  indR = duplicated(dtemp[(floor(M/2)+1):M])
  dtemp = dtemp[!c(indL, indR)]
  dens_temp = exp(-lqd[!c(indL, indR)]);
  
  # Interpolate to dSup and normalize
  dSup = seq(dtemp[1], dtemp[length(dtemp)], length.out = M)
  dens = approx(x = dtemp, y = dens_temp, xout = dSup, rule = c(2,2))[[2]]
  dens = dens/trapzRcpp(X = dSup,Y = dens)*(lqdSup[M] - lqdSup[1]); # Normalize, accounting for boundary cutoff
  
  return(list('dSup' = dSup, 'dens' = dens))
}

Cal.Summary <- function(Est = NA, Pval = NA){
  
  if (!anyNA(Est)){
    iSummary <- c(mean(Est), std(Est), quantile(Est,c(0.25,0.5,0.75)))
  } else{
    iSummary <- NULL
  }
  
  if (!anyNA(Pval)){
    rej1 <- mean(Pval<0.01)
    rej5 <- mean(Pval<0.05)
    rej10 <- mean(Pval<0.10)
    
    iSummary <- c(iSummary, rej1, rej5, rej10)
  }
  
  iSummary
  
}

My.CreateScreePlot <- function(fpcaObj, ...){ 
  
  args1 <- list( main="Scree plot", ylab='% of variance explained', xlab='Number of components')  
  inargs <- list(...)
  args1[names(inargs)] <- inargs
  
  ys <- fpcaObj$cumFVE ;
  
  
  if( !is.vector(ys) ){ 
    stop('Please use a vector as input.')   
  }
  if(max(ys) > 1){
    warning('The maximum number in the input vector is larger than 1; are you sure it is right?');
  }
  if(any(ys < 0) || any(diff(ys) <0) ){
    stop('This is not a valid cumulative FVE vector.')
  }
  
  dfbar <- do.call( barplot, c( args1, list( ylim=c(0,1.05)), list(axes=FALSE), list(height =  rep(NA,length(ys))) ) )
  
  abline(h=(seq(0,1,.05)), col="lightgray", lty="dotted")
  barplot(c(ys[1], diff(ys)), add = TRUE , names.arg = as.character(1:fpcaObj$selectK))
  lines(dfbar, y= ys, col='red')
  points(dfbar, y= ys, col='red')
  legend("right", "CPV", col='red', lty=1, pch=1, bty='n') 
  
}

Make.Reg.Data.h <- function(CPI, UNEMP, Scores, h){
  
  hplus0     <- paste0('plus', h)
  hplus1     <- paste0('plus', h + 1)
  
  myData.Temp <- cbind(CPI, UNEMP[[hplus0]]) %>% rename( UNEMP.hplus0 = "UNEMP[[hplus0]]")
  scores.lag <- dplyr::lag(Scores[[hplus1]], n=1)
  PCS <- data.frame(UniDate, scores.lag)
  colnames(PCS) <- c("Date", paste0("score", 1:3))
  myRegData.h <-  merge(myData.Temp, PCS, by.x="Date", by.y="Date") %>% 
    dplyr::select(Date, ID, contains(hplus0), contains(hplus1), UNEMP.hplus0, score1, score2, score3) %>% 
    rename( CPI.hplus0 = all_of(hplus0), CPI.hplus1 = all_of(hplus1))
  
  myRegData.h
}

Make.Reg.Data.Moments <- function(CPI, UNEMP, Scores, Moments, h){
  
  hplus0     <- paste0('plus', h)
  hplus1     <- paste0('plus', h + 1)
  
  myData.Temp <- cbind(CPI, UNEMP[[hplus0]]) %>% rename( UNEMP.hplus0 = "UNEMP[[hplus0]]")
  scores.lag  <- dplyr::lag(Scores[[hplus1]], n=1)
  Moments.lag <- dplyr::lag(Moments[[hplus1]], n=1)
  
  PCS <- data.frame(UniDate, scores.lag, Moments.lag)
  colnames(PCS) <- c("Date", paste0("score", 1:3), "Mean", "Median","SD","IQR","Skew","Kurtosis")
  myRegData.h <-  merge(myData.Temp, PCS, by.x="Date", by.y="Date") %>% 
    dplyr::select(Date, ID, contains(hplus0), contains(hplus1), UNEMP.hplus0, score1, score2, score3, Mean, Median, SD, IQR, Skew, Kurtosis) %>% 
    rename( CPI.hplus0 = all_of(hplus0), CPI.hplus1 = all_of(hplus1))
  
  myRegData.h
}

Make.Reg.Data.External <- function(CPI, UNEMP, Scores, External, h){
  
  hplus0     <- paste0('plus', h)
  hplus1     <- paste0('plus', h + 1)
  
  myData.Temp <- cbind(CPI, UNEMP[[hplus0]]) %>% rename( UNEMP.hplus0 = "UNEMP[[hplus0]]")
  scores.lag  <- dplyr::lag(Scores[[hplus1]], n=1)
  
  PCS <- data.frame(UniDate, scores.lag, External$VIX, External$Recession, External$Inflation)
  colnames(PCS) <- c("Date", paste0("score", 1:3), "VIX", "Recession", "Inflation")
  myRegData.h <-  merge(myData.Temp, PCS, by.x="Date", by.y="Date") %>% 
    dplyr::select(Date, ID, contains(hplus0), contains(hplus1), UNEMP.hplus0, score1, score2, score3, VIX, Recession, Inflation) %>% 
    rename( CPI.hplus0 = all_of(hplus0), CPI.hplus1 = all_of(hplus1))
  
  myRegData.h
}

Make.Reg.Data.Moments.External <- function(CPI, UNEMP, Scores, Moments, External, h){
  
  hplus0     <- paste0('plus', h)
  hplus1     <- paste0('plus', h + 1)
  
  myData.Temp <- cbind(CPI, UNEMP[[hplus0]]) %>% rename( UNEMP.hplus0 = "UNEMP[[hplus0]]")
  scores.lag  <- dplyr::lag(Scores[[hplus1]], n=1)
  Moments.lag <- dplyr::lag(Moments[[hplus1]], n=1)
  
  PCS <- data.frame(UniDate, scores.lag, Moments.lag, External$VIX, External$Recession, External$Inflation)
  colnames(PCS) <- c("Date", paste0("score", 1:3), "Mean", "Median","SD","IQR","Skew","Kurtosis", "VIX", "Recession", "Inflation")
  myRegData.h <-  merge(myData.Temp, PCS, by.x="Date", by.y="Date") %>% 
    dplyr::select(Date, ID, contains(hplus0), contains(hplus1), UNEMP.hplus0, score1, score2, score3, Mean, Median, SD, IQR, Skew, Kurtosis, VIX, Recession, Inflation) %>% 
    rename( CPI.hplus0 = all_of(hplus0), CPI.hplus1 = all_of(hplus1))
  
  myRegData.h
}

Make.Reg.Data.hybrid <- function(CPI, UNEMP, Scores, h){
  
  hplus0     <- paste0('plus', h)
  hplus1     <- paste0('plus', h + 1)
  hminus1    <- paste0('plus', h - 1)
  
  myData.Temp <- cbind(CPI, UNEMP[[hplus0]]) %>% rename( UNEMP.hplus0 = "UNEMP[[hplus0]]")
  scores.lag <- dplyr::lag(Scores[[hplus1]], n=1)
  PCS <- data.frame(UniDate, scores.lag)
  colnames(PCS) <- c("Date", paste0("score", 1:3))
  myRegData.h <-  merge(myData.Temp, PCS, by.x="Date", by.y="Date") %>% 
    dplyr::select(Date, ID, contains(hplus0), contains(hplus1), contains(hminus1), UNEMP.hplus0, score1, score2, score3) %>% 
    rename( CPI.hplus0 = all_of(hplus0), CPI.hplus1 = all_of(hplus1),  CPI.hminus1 = all_of(hminus1))
  
  myRegData.h
}

Make.Reg.Data.hybrid.2Lag <- function(CPI, UNEMP, Scores, h){
  
  hplus0     <- paste0('plus', h)
  hplus1     <- paste0('plus', h + 1)
  hminus1    <- paste0('plus', h - 1)
  hminus2    <- paste0('plus', h - 2)
  
  myData.Temp <- cbind(CPI, UNEMP[[hplus0]]) %>% rename( UNEMP.hplus0 = "UNEMP[[hplus0]]")
  scores.lag <- dplyr::lag(Scores[[hplus1]], n=1)
  PCS <- data.frame(UniDate, scores.lag)
  colnames(PCS) <- c("Date", paste0("score", 1:3))
  myRegData.h <-  merge(myData.Temp, PCS, by.x="Date", by.y="Date") %>% 
    dplyr::select(Date, ID, contains(hplus0), contains(hplus1), contains(hminus1), contains(hminus2), UNEMP.hplus0, score1, score2, score3) %>% 
    rename( CPI.hplus0 = all_of(hplus0), CPI.hplus1 = all_of(hplus1),  CPI.hminus1 = all_of(hminus1), CPI.hminus2 = all_of(hminus2))
  
  myRegData.h
}

Make.Reg.Data.decompseS2 <- function(CPI, UNEMP, Scores, h){
  
  hplus0     <- paste0('plus', h)
  hplus1     <- paste0('plus', h + 1)
  
  myData.Temp <- cbind(CPI, UNEMP[[hplus0]]) %>% rename( UNEMP.hplus0 = "UNEMP[[hplus0]]")
  scores.lag <- dplyr::lag(Scores[[hplus1]], n=1)
  PCS <- data.frame(UniDate, scores.lag)
  colnames(PCS) <- c("Date", paste0("score", 1:3))
  myRegData.h <-  merge(myData.Temp, PCS, by.x="Date", by.y="Date") %>% 
    dplyr::select(Date, ID, contains(hplus0), contains(hplus1), UNEMP.hplus0, score1, score2, score3) %>% 
    rename( CPI.hplus0 = all_of(hplus0), CPI.hplus1 = all_of(hplus1))
  
  myRegData.h$score2_pos <- pmax(myRegData.h$score2,0)
  myRegData.h$score2_neg <- pmin(0,myRegData.h$score2)
  myRegData.h
}

Make.Reg.Data.Seperately <- function(CPI, UNEMP, scores.hplus1, h){
  
  hplus0 <- paste0('plus', h)
  hplus1 <- paste0('plus', h + 1)
  
  myData.Temp <- cbind(CPI, UNEMP[[hplus0]]) %>% rename( UNEMP.hplus0 = "UNEMP[[hplus0]]")
  scores.lag  <- dplyr::lag(scores.hplus1, n=1)
  PCS <- data.frame(UniDate, scores.lag)
  colnames(PCS) <- c("Date", paste0("score", 1:3))
  myRegData.h <-  merge(myData.Temp, PCS, by.x="Date", by.y="Date") %>% 
    dplyr::select(Date, ID, contains(hplus0), contains(hplus1), UNEMP.hplus0, score1, score2, score3) %>% 
    rename( CPI.hplus0 = all_of(hplus0), CPI.hplus1 = all_of(hplus1))
  
  myRegData.h
}

Cal.Charc <- function(CPI, Res.ind){
  Resp <- Res.ind[,'ID']
  nResp <- length(Resp)
  Charc <- matrix(NA, nrow = nResp, ncol = 3)
  
  for (i in 1:nResp){
    iFore <- CPI %>% subset(ID == Resp[i]) %>% dplyr::select("plus0","plus1","plus2","plus3","Date")
    
    iCut <- as.yearqtr(1981.50 + 73/4)
    iEarly <-  subset(iFore, Date <= iCut)
    iLate <- subset(iFore, iCut < Date)
    iEarly.nFore <- sum(apply(iEarly, 1, function(x) all(!is.na(x))))
    iLate.nFore <- sum(apply(iLate, 1, function(x) all(!is.na(x))))
    iPeriod <- iEarly.nFore/(iEarly.nFore + iLate.nFore)
    
    iFore  <- dplyr::select(iFore, -"Date")
    iSigma <- mean(apply(iFore, 2, sd, na.rm = TRUE)) # sd(iFore, na.rm = TRUE)
    
    Charc[i, 1] <- iPeriod
    Charc[i, 2] <- iSigma
    Charc[i, 3] <- sign(Res.ind[i,"UNEMP"])
  }
  
  colnames(Charc) <- c("Period_j","Sigma_tilde","Sign_UNEMP")
  
  Charc <- cbind(Res.ind, Charc)
  data.frame(Charc)
}

Count.Pairs <- function(JTest){
  JTest.Summary <- matrix(NA, 4,3)
  Alpha <- c(0.01, 0.05, 0.10)
  nAlpha <- length(Alpha)
  
  for (i in 1:nAlpha){
    JTest.Summary[1,i] <- mean(apply((JTest < Alpha[i]), 1, identical, c(TRUE,FALSE)))
    JTest.Summary[2,i] <- mean(apply((JTest < Alpha[i]), 1, identical, c(FALSE,TRUE)))
    JTest.Summary[3,i] <- mean(apply((JTest < Alpha[i]), 1, identical, c(TRUE,TRUE)))
    JTest.Summary[4,i] <- mean(apply((JTest < Alpha[i]), 1, identical, c(FALSE,FALSE)))
  }

  row.names(JTest.Summary) <- c("FPCR encompasses the moments-based model","The moments-based model encompasses FPCR","Both models encompass each other", "Both models do not encompass each other")
  colnames(JTest.Summary) <- c("rej_1", "rej_5", "rej_10")
  JTest.Summary
}







