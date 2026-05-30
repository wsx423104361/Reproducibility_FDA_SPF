# Preliminary--------------------------------------------------------------
rm(list = ls())
Sys.setlocale(category = "LC_ALL", locale = "English")
par(mfrow = c(1, 1))

# Library -----------------------------------------------------------------
library(readxl)
library(zoo)
library(pracma)
library(e1071)
library(fdadensity)
library(fdapace)
library(plm)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(ftsa)
library(lmtest)
library(tidyverse)
source("myFunctions.R")

# Read and Clean Data -----------------------------------------------------
CPI <- read_excel("./Data/Individual_CPI_2022Q2.xlsx", col_types = "numeric")
UNEMP <- read_excel("./Data/Individual_UNEMP_2022Q2.xlsx", col_types = "numeric")
CPI.long <- read_excel("./Data/CPI10_BlueChip.xlsx", col_types = "numeric")

CPI.Mean <- read_excel("./Data/Aggregate_CPI_2022Q2.xlsx", sheet = "Mean")
CPI.Median <- read_excel("./Data/Aggregate_CPI_2022Q2.xlsx", sheet = "Median")
CPI.D1 <- read_excel("./Data/Aggregate_CPI_2022Q2.xlsx", sheet = "Dispersion1")

External <- read_excel("./Data/External_Variables.xlsx")

## Time Range between 1981Q3 to 2022Q2
CPI$Date <- as.yearqtr(CPI$YEAR + (CPI$QUARTER-1)/4)
UNEMP$Date <- as.yearqtr(UNEMP$YEAR + (UNEMP$QUARTER-1)/4)
CPI.long$Date <- as.yearqtr(CPI.long$YEAR + (CPI.long$QUARTER-1)/4)

UniDate <- as.yearqtr(1981.50 + seq(0, 163)/4)
RawDate <- CPI$Date
CPI.long <- subset(CPI.long, Date %in% UniDate)

CPI <- subset(CPI, Date %in% UniDate)
UNEMP <- subset(UNEMP, Date %in% UniDate)

Hname <- paste0('plus', 0:4)
nH <- length(Hname)

# Obtain Density Curves for Figure 1 ------------------------------------------------
nDate <- length(UniDate)
Density <- list()
npoint0 <- 512
for (t in 1:nDate){
  qtr <- UniDate[t]
  CPI.qtr <- subset(CPI, Date %in% qtr)
  Density$plus1[[t]] <- density(CPI.qtr$plus1, bw = "nrd0", from=-5, to=10, n = npoint0, na.rm=TRUE) %>% normalize.density # Note: this setting is for a better visualization of the density curves
}

Density$plus1 <- matrix(unlist(Density$plus1),npoint0,nDate) 
write.table(Density$plus1, "Data/Density_h1.csv", sep = ",", col.names = FALSE, row.names = FALSE)

# Obtain LQD Curves for Figure 1 ------------------------------------------------
npoint <- 512/2
bw.choice <- "nrd0" # SJ; nrd0
lqdSup <- seq(0, 1, length.out = npoint)
alpha <- 0.01 

Density <- list()
for (t in 1:nDate){
  t
  qtr <- UniDate[t]
  CPI.qtr <- subset(CPI, Date %in% qtr)
  
  for (h in 1:nH){
    h
    hname <- Hname[h]
    CPI.qtr.h <- CPI.qtr[[hname]]
    
    t.h.dens <- density(CPI.qtr.h, bw = bw.choice, from=min(CPI.qtr.h, na.rm = T), to=max(CPI.qtr.h, na.rm = T), n = npoint, na.rm=TRUE)
    t.h.from <- quantile(CPI.qtr.h, 0.01, na.rm = T) # 0.01
    t.h.to <- quantile(CPI.qtr.h, 0.99, na.rm = T) # 0.99
    t.h.dens <- density(CPI.qtr.h, bw = bw.choice, from=t.h.from, to=t.h.to, n = npoint, na.rm=TRUE)
    t.h.dens.dSup <- t.h.dens$x
    t.h.dens.norm <- t.h.dens %>% normalize.density() %>% RegulariseByAlpha(x=t.h.dens$x, alpha=alpha)
    
    Density[[hname]][['dens']][[t]] <- t.h.dens.norm
    Density[[hname]][['qd']][[t]]   <- dens2qd(dens=t.h.dens.norm, dSup = t.h.dens.dSup, qdSup =lqdSup)
    Density[[hname]][['lqd']][[t]]  <- dens2lqd(dens=t.h.dens.norm,dSup = t.h.dens.dSup, lqdSup=lqdSup)
    Density[[hname]][['dSup']][[t]] <- t.h.dens.dSup
    
    Density[[hname]][['Mean']][[t]] <- mean(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['Median']][[t]] <- median(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['SD']][[t]] <- sd(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['IQR']][[t]] <- IQR(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['Skew']][[t]] <- skewness(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['Kurt']][[t]] <- kurtosis(CPI.qtr.h, na.rm=TRUE) + 3
    Density[[hname]][['Curv']][[t]] <- trapz(lqdSup, c(0,0,diff(diff(Density[[hname]][['lqd']][[t]])))^2)
  }
}
write.table(Density[['plus1']]$lqd, "Data/LQD_h1.csv", sep = ",", col.names = FALSE, row.names = FALSE)

# Table 1 Descriptive Statistics -------------------------------------
DS.CPI <- cal.DS(CPI,nDate,UniDate)
DS.UNEMP <- cal.DS(UNEMP,nDate,UniDate)
write.csv(DS.CPI$DS.summary, "Tables/Table01_Summ_Stat_CPI_Upper_Panel.csv")
write.csv(DS.UNEMP$DS.summary, "Tables/Table01_Summ_Stat_Unemployment_Lower_Panel.csv")

# FPCA by lqd ---------------------------------------------------------
lqd.plus1 <- matrix(unlist(Density[['plus1']]$lqd), nrow=nDate, ncol=npoint, byrow = T)
lqd.plus2 <- matrix(unlist(Density[['plus2']]$lqd), nrow=nDate, ncol=npoint, byrow = T)
lqd.plus3 <- matrix(unlist(Density[['plus3']]$lqd), nrow=nDate, ncol=npoint, byrow = T)
lqd.plus4 <- matrix(unlist(Density[['plus4']]$lqd), nrow=nDate, ncol=npoint, byrow = T)
lqd.pool  <- rbind(lqd.plus1, lqd.plus2, lqd.plus3, lqd.plus4)

LQD <- MakeFPCAInputs(tVec = lqdSup, yVec = lqd.pool)
FPCA.lqd <- fdapace::FPCA(Ly = LQD$Ly, Lt = LQD$Lt)

# Figure 2
setEPS()
postscript("Plots/Fig02_No_Components.eps", width = 8, height = 6)
My.CreateScreePlot(FPCA.lqd, ylab = "% of variance explained", main = "")
dev.off()

# Extract scores
Mu <- FPCA.lqd$mu
PC <- FPCA.lqd$phi[,1:3]
scores.allh <- FPCA.lqd$xiEst[,1:3]

# Figure 3
setEPS()
postscript("Plots/Fig03_EFPCs.eps", width = 8, height = 6)
matplot(FPCA.lqd$workGrid, PC, type='l', xlab="", ylab="EFPC", lwd = 2, col = c("black","red","blue"), lty = 1:3); grid(); legend("top",c("EFPC 1","EFPC 2","EFPC 3"), cex=.8, col = c("black","red","blue"), lty = 1:3)
dev.off()

# Figure 4
setEPS()
postscript("Plots/Fig04_Scores.eps", width = 20, height = 15)
par(mfrow=c(2,2))
Scores <- list()
Scores$plus1 <- scores.allh[(nDate*0+1):(nDate*1),]
Scores$plus2 <- scores.allh[(nDate*1+1):(nDate*2),]
Scores$plus3 <- scores.allh[(nDate*2+1):(nDate*3),]
Scores$plus4 <- scores.allh[(nDate*3+1):(nDate*4),]
matplot(UniDate, Scores$plus1, type='l', xlab="Date", ylab="scores", main = "h=0", lwd = 1.5, col = c("black","red","blue"), lty = 1:3); grid(); legend("top",c("score 1","score 2","score 3"), col = c("black","red","blue"), lty = 1:3)
matplot(UniDate, Scores$plus2, type='l', xlab="Date", ylab="scores", main = "h=1", lwd = 1.5, col = c("black","red","blue"), lty = 1:3); grid(); legend("top",c("score 1","score 2","score 3"), col = c("black","red","blue"), lty = 1:3)
matplot(UniDate, Scores$plus3, type='l', xlab="Date", ylab="scores", main = "h=2", lwd = 1.5, col = c("black","red","blue"), lty = 1:3); grid(); legend("top",c("score 1","score 2","score 3"), col = c("black","red","blue"), lty = 1:3)
matplot(UniDate, Scores$plus4, type='l', xlab="Date", ylab="scores", main = "h=3", lwd = 1.5, col = c("black","red","blue"), lty = 1:3); grid(); legend("top",c("score 1","score 2","score 3"), col = c("black","red","blue"), lty = 1:3)
dev.off()

# Verifying the Principle Component Scores
Mean.plus1 <- unlist(Density[['plus1']][['Mean']])
Mean.plus2 <- unlist(Density[['plus2']][['Mean']])
Mean.plus3 <- unlist(Density[['plus3']][['Mean']])
Mean.plus4 <- unlist(Density[['plus4']][['Mean']])
Mean.pool <- c(Mean.plus1, Mean.plus2, Mean.plus3, Mean.plus4)

Median.plus1 <- unlist(Density[['plus1']][['Median']])
Median.plus2 <- unlist(Density[['plus2']][['Median']])
Median.plus3 <- unlist(Density[['plus3']][['Median']])
Median.plus4 <- unlist(Density[['plus4']][['Median']])
Median.pool <- c(Median.plus1, Median.plus2, Median.plus3, Median.plus4)

SD.plus1 <- unlist(Density[['plus1']][['SD']])
SD.plus2 <- unlist(Density[['plus2']][['SD']])
SD.plus3 <- unlist(Density[['plus3']][['SD']])
SD.plus4 <- unlist(Density[['plus4']][['SD']])
SD.pool <- c(SD.plus1, SD.plus2, SD.plus3, SD.plus4)

IQR.plus1 <- unlist(Density[['plus1']][['IQR']])
IQR.plus2 <- unlist(Density[['plus2']][['IQR']])
IQR.plus3 <- unlist(Density[['plus3']][['IQR']])
IQR.plus4 <- unlist(Density[['plus4']][['IQR']])
IQR.pool <- c(IQR.plus1, IQR.plus2, IQR.plus3, IQR.plus4)

Skew.plus1 <- unlist(Density[['plus1']][['Skew']])
Skew.plus2 <- unlist(Density[['plus2']][['Skew']])
Skew.plus3 <- unlist(Density[['plus3']][['Skew']])
Skew.plus4 <- unlist(Density[['plus4']][['Skew']])
Skew.pool <- c(Skew.plus1, Skew.plus2, Skew.plus3, Skew.plus4)

Kurt.plus1 <- unlist(Density[['plus1']][['Kurt']])
Kurt.plus2 <- unlist(Density[['plus2']][['Kurt']])
Kurt.plus3 <- unlist(Density[['plus3']][['Kurt']])
Kurt.plus4 <- unlist(Density[['plus4']][['Kurt']])
Kurt.pool <- c(Kurt.plus1, Kurt.plus2, Kurt.plus3, Kurt.plus4)

Curv.plus1 <- unlist(Density[['plus1']][['Curv']])
Curv.plus2 <- unlist(Density[['plus2']][['Curv']])
Curv.plus3 <- unlist(Density[['plus3']][['Curv']])
Curv.plus4 <- unlist(Density[['plus4']][['Curv']])
Curv.pool <- c(Curv.plus1, Curv.plus2, Curv.plus3, Curv.plus4)

Oil.pool <- rep(External$Oil, times = 4)
Lab.pool <- rep(External$LabourCost, times = 4)

# Table 2
EFPC.Cor <- matrix(nrow = 8, ncol = 3)
EFPC.Cor[1,] <- cor(scores.allh, Mean.pool)
EFPC.Cor[2,] <- cor(scores.allh, Median.pool)
EFPC.Cor[3,] <- cor(scores.allh, SD.pool)
EFPC.Cor[4,] <- cor(scores.allh, IQR.pool)
EFPC.Cor[5,] <- cor(scores.allh, Skew.pool)
EFPC.Cor[6,] <- cor(scores.allh, Kurt.pool)
EFPC.Cor[7,] <- cor(scores.allh, Oil.pool)
EFPC.Cor[8,] <- cor(scores.allh, Lab.pool)
row.names(EFPC.Cor) <- c("Mean", "Median","SD", "IQR", "Skewness", "Kurtosis", "Oil_Price", "Labour_Cost")
colnames(EFPC.Cor) <- c("Score1","Score2","Score3")
write.csv(EFPC.Cor, "Tables/Table02_Corr_Pearson_Left_Panel.csv")

EFPC.rho <- matrix(nrow = 8, ncol = 3)
EFPC.rho[1,] <- cor(scores.allh, Mean.pool, method = "spearman")
EFPC.rho[2,] <- cor(scores.allh, Median.pool, method = "spearman")
EFPC.rho[3,] <- cor(scores.allh, SD.pool, method = "spearman")
EFPC.rho[4,] <- cor(scores.allh, IQR.pool, method = "spearman")
EFPC.rho[5,] <- cor(scores.allh, Skew.pool, method = "spearman")
EFPC.rho[6,] <- cor(scores.allh, Kurt.pool, method = "spearman")
EFPC.rho[7,] <- cor(scores.allh, Oil.pool, method = "spearman")
EFPC.rho[8,] <- cor(scores.allh, Lab.pool, method = "spearman")
row.names(EFPC.rho) <- c("Mean", "Median","SD", "IQR", "Skewness", "Kurtosis", "Oil_Price", "Labour_Cost")
colnames(EFPC.rho) <- c("Score1","Score2","Score3")
write.csv(EFPC.rho, "Tables/Table02_Spearman_right_Panel.csv")

# Figure 5
cairo_pdf("Plots/Fig05_Scores_Moments.pdf", width = 18, height = 6)
par(mfrow=c(1,3), xpd = FALSE, mar = c(6, 6, 2, 1) + 0.1, family = "serif")
plot(SD.pool, scores.allh[,1],   ylab="First Score",  xlab = "SD", pch=18, cex.lab = 1.7, cex.axis = 2); abline(lm(scores.allh[,1]~SD.pool), col="red"); lines(lowess(SD.pool,scores.allh[,1]), col="blue"); grid()
plot(Skew.pool, scores.allh[,2],  ylab="Second Score", xlab = "skewness", pch=18, cex.lab = 1.7, cex.axis = 2); abline(lm(scores.allh[,2]~Skew.pool), col="red"); lines(lowess(Skew.pool,scores.allh[,2]), col="blue"); grid()
plot(Kurt.pool, scores.allh[,3], ylab="Third Score",  xlab = "kurtosis", pch=18, cex.lab = 1.7, cex.axis = 2); abline(lm(scores.allh[,3]~Kurt.pool), col="red"); lines(lowess(Kurt.pool,scores.allh[,3]), col="blue"); grid()
dev.off()

# Baseline of Functional Linear Model  ------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=0)
myRegData.plus1 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=1)
myRegData.plus2 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=2)
myRegData.plus3 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

# Individual Regression
uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
Res2 <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    iRes.reduce.summary <- summary(iRes.reduce)
    
    Res[[i]]  <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    Res2[[i]] <- c(iID, matrix(t(iRes.reduce.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.reduce.summary$adj.r.squared, nobs(iRes.reduce))
    
    no.col <- length(Res[[i]])
    no.col2 <- length(Res2[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

# Figure 6
Fig06 <- Res.DF %>%
  dplyr::select(c("CPI","score1","score2","score3","UNEMP")) %>% gather() %>%
  ggplot( aes(x=key, y=value, fill=key)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9, position = position_jitter(seed = 123)) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + ggtitle(" ") + xlab(" ") + ylab(" ")
ggsave(filename = "Plots/Fig06_Boxplot.pdf", plot = Fig06, width = 14, height = 8, device = cairo_pdf)

# Table 3 (Upper Panel)
Res.ind.Summary <- matrix(nrow=9,ncol=8)
Res.ind.Summary[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary) <- c("Intercept","CPI","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary, "Tables/Table03_Baseline_Upper_Panel.csv")

# Regression without Scores 1-3 
myFormula2 <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

Res2 <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula2, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    Res2[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res2[[i]])
  }
}

# Table 3 (Lower Panel)
Res2.ind <- matrix(unlist(Res2), ncol=no.col2, byrow=T)
colnames(Res2.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","UNEMP","UNEMP_pval","Adj_R2", "No_Obs")
Res2.ind.Summary <- matrix(nrow=5,ncol=8)
Res2.ind.Summary[1,] <- Cal.Summary(Res2.ind[,"Intercept"], Res2.ind[,"Intercept_pval"])
Res2.ind.Summary[2,] <- Cal.Summary(Res2.ind[,"CPI"], Res2.ind[,"CPI_pval"])
Res2.ind.Summary[3,] <- Cal.Summary(Res2.ind[,"UNEMP"], Res2.ind[,"UNEMP_pval"])
Res2.ind.Summary[4,1:5] <- Cal.Summary(Res2.ind[,"Adj_R2"])
Res2.ind.Summary[5,1:5] <- Cal.Summary(Res2.ind[,"No_Obs"])
row.names(Res2.ind.Summary) <- c("Intercept","CPI", "UNEMP","Adj_R2", "No_Obs" )
colnames(Res2.ind.Summary) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res2.ind.Summary, "Tables/Table03_Baseline_Lower_Panel.csv")

# Characteristics of Respondents ------------------------------------------
Charc <- Cal.Charc(CPI, Res.ind)
Crit <- c("score1_pval","score2_pval","score3_pval","F_pval")
nCrit <- length(Crit)
Res.Charc <- list() 

for (i in 1:nCrit) {
  iRes.Charc <- matrix(NA, nrow = 7, ncol = 5)
  
  iCrit <- Crit[i]
  Sig.idx <- Charc[[iCrit]] < 0.05
  
  Charc.sel <- dplyr::select(Charc,c("UNEMP", "Sign_UNEMP", "CPI","Period_j","Sigma_tilde"))
  Sig.Group <- Charc.sel[Sig.idx, ]
  InSig.Group <- Charc.sel[!Sig.idx, ]
  
  iRes.Charc[1, ] <- apply(Sig.Group, 2, mean)
  iRes.Charc[2, ] <- apply(Sig.Group, 2, sd)  
  iRes.Charc[3, ] <- apply(Sig.Group, 2, length) 
  iRes.Charc[4, ] <- apply(InSig.Group, 2, mean)
  iRes.Charc[5, ] <- apply(InSig.Group, 2, sd)
  iRes.Charc[6, ] <- apply(InSig.Group, 2, length)
  
  iRes.Charc[7, 1] <- t.test(Sig.Group$UNEMP, InSig.Group$UNEMP)$p.value
  iRes.Charc[7, 2] <- t.test(Sig.Group$Sign_UNEMP, InSig.Group$Sign_UNEMP)$p.value
  iRes.Charc[7, 3] <- t.test(Sig.Group$CPI, InSig.Group$CPI)$p.value
  iRes.Charc[7, 4] <- t.test(Sig.Group$Period_j, InSig.Group$Period_j)$p.value
  iRes.Charc[7, 5] <- t.test(Sig.Group$Sigma_tilde, InSig.Group$Sigma_tilde)$p.value
  
  Res.Charc[[i]] <- iRes.Charc
}

Res.Charc <- do.call(rbind, Res.Charc)
colnames(Res.Charc) <- c("alpha_j","sign_alpha_j","beta_j","Period_j","sigma_j")
rownames(Res.Charc) <- c("Score1_Sig_Mean","Score1_Sig_sd","Score1_Sig_No","Score1_Insig_Mean","Score1_InSig_sd","Score1_InSig_No","Score1_t_test_pval",
                         "Score2_Sig_Mean","Score2_Sig_sd","Score2_Sig_No","Score2_Insig_Mean","Score2_InSig_sd","Score2_InSig_No","Score2_t_test_pval",
                         "Score3_Sig_Mean","Score3_Sig_sd","Score3_Sig_No","Score3_Insig_Mean","Score3_InSig_sd","Score3_InSig_No","Score3_t_test_pval",
                         "F_test_Sig_Mean","F_test_Sig_sd","F_test_Sig_No","F_test_Insig_Mean","F_test_InSig_sd","F_test_InSig_No","F_test_t_test_pval")
write.csv(Res.Charc, "Tables/Table04_Characteristics.csv")

# Include Mean ----------------------------------------------------------------------------------------------------------------------------------
Moments <- list()
Moments$plus1 <- cbind(Mean.plus1, Median.plus1, SD.plus1, IQR.plus1, Skew.plus1, Kurt.plus1)
Moments$plus2 <- cbind(Mean.plus2, Median.plus2, SD.plus2, IQR.plus2, Skew.plus2, Kurt.plus2)
Moments$plus3 <- cbind(Mean.plus3, Median.plus3, SD.plus3, IQR.plus3, Skew.plus3, Kurt.plus3)
Moments$plus4 <- cbind(Mean.plus4, Median.plus4, SD.plus4, IQR.plus4, Skew.plus4, Kurt.plus4)

myRegData.plus0 <- Make.Reg.Data.Moments(CPI, UNEMP, Scores, Moments, h=0)
myRegData.plus1 <- Make.Reg.Data.Moments(CPI, UNEMP, Scores, Moments, h=1)
myRegData.plus2 <- Make.Reg.Data.Moments(CPI, UNEMP, Scores, Moments, h=2)
myRegData.plus3 <- Make.Reg.Data.Moments(CPI, UNEMP, Scores, Moments, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + Mean + score1 + score2 + score3 + UNEMP.hplus0') 
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1  + Mean + UNEMP.hplus0') 

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
Res2 <- list()

for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(Mean)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    iRes.reduce.summary <- summary(iRes.reduce)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    Res2[[i]] <- c(iID, matrix(t(iRes.reduce.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.reduce.summary$adj.r.squared, nobs(iRes.reduce))
    
    no.col <- length(Res[[i]])
    no.col2 <- length(Res2[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval", "Mean", "Mean_pval", "score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","Mean_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.IncMean <- matrix(nrow=10,ncol=8)
Res.ind.Summary.IncMean[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.IncMean[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.IncMean[3,] <- Cal.Summary(Res.ind[,"Mean"], Res.ind[,"Mean_pval"])
Res.ind.Summary.IncMean[4,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.IncMean[5,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.IncMean[6,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.IncMean[7,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.IncMean[8,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.IncMean[9,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.IncMean[10,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.IncMean) <- c("Intercept","CPI", "Mean","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.IncMean) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.IncMean, "Tables/Table05_Include_Mean.csv")

# Include Higher Moments-------------------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.Moments(CPI, UNEMP, Scores, Moments, h=0)
myRegData.plus1 <- Make.Reg.Data.Moments(CPI, UNEMP, Scores, Moments, h=1)
myRegData.plus2 <- Make.Reg.Data.Moments(CPI, UNEMP, Scores, Moments, h=2)
myRegData.plus3 <- Make.Reg.Data.Moments(CPI, UNEMP, Scores, Moments, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

# Include SD
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + SD + score1 + score2 + score3 + UNEMP.hplus0') 
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1  + SD + UNEMP.hplus0') 

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
Res2 <- list()

for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(Mean)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    iRes.reduce.summary <- summary(iRes.reduce)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    Res2[[i]] <- c(iID, matrix(t(iRes.reduce.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.reduce.summary$adj.r.squared, nobs(iRes.reduce))
    
    no.col <- length(Res[[i]])
    no.col2 <- length(Res2[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval", "SD", "SD_pval", "score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","SD_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.IncSD <- matrix(nrow=10,ncol=8)
Res.ind.Summary.IncSD[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.IncSD[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.IncSD[3,] <- Cal.Summary(Res.ind[,"SD"], Res.ind[,"SD_pval"])
Res.ind.Summary.IncSD[4,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.IncSD[5,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.IncSD[6,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.IncSD[7,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.IncSD[8,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.IncSD[9,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.IncSD[10,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.IncSD) <- c("Intercept","CPI", "SD","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.IncSD) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.IncSD, "Tables/Table06_Include_SD_Upper_Panel.csv")

# Include Skewness
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + Skew + score1 + score2 + score3 + UNEMP.hplus0') 
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1  + Skew + UNEMP.hplus0') 

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
Res2 <- list()

for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(Mean)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    iRes.reduce.summary <- summary(iRes.reduce)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    Res2[[i]] <- c(iID, matrix(t(iRes.reduce.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.reduce.summary$adj.r.squared, nobs(iRes.reduce))
    
    no.col <- length(Res[[i]])
    no.col2 <- length(Res2[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval", "Skew", "Skew_pval", "score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","Skew_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.IncSkew <- matrix(nrow=10,ncol=8)
Res.ind.Summary.IncSkew[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.IncSkew[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.IncSkew[3,] <- Cal.Summary(Res.ind[,"Skew"], Res.ind[,"Skew_pval"])
Res.ind.Summary.IncSkew[4,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.IncSkew[5,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.IncSkew[6,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.IncSkew[7,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.IncSkew[8,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.IncSkew[9,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.IncSkew[10,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.IncSkew) <- c("Intercept","CPI", "Skew","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.IncSkew) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.IncSkew, "Tables/Table06_Include_Skew_Middle_Panel.csv")


# Include Kurtosis
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + Kurtosis + score1 + score2 + score3 + UNEMP.hplus0') 
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1  + Kurtosis + UNEMP.hplus0') 

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
Res2 <- list()

for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(Mean)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    iRes.reduce.summary <- summary(iRes.reduce)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    Res2[[i]] <- c(iID, matrix(t(iRes.reduce.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.reduce.summary$adj.r.squared, nobs(iRes.reduce))
    
    no.col <- length(Res[[i]])
    no.col2 <- length(Res2[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval", "Kurt", "Kurt_pval", "score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","Kurt_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.IncKurt <- matrix(nrow=10,ncol=8)
Res.ind.Summary.IncKurt[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.IncKurt[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.IncKurt[3,] <- Cal.Summary(Res.ind[,"Kurt"], Res.ind[,"Kurt_pval"])
Res.ind.Summary.IncKurt[4,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.IncKurt[5,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.IncKurt[6,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.IncKurt[7,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.IncKurt[8,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.IncKurt[9,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.IncKurt[10,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.IncKurt) <- c("Intercept","CPI", "Kurt","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.IncKurt) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.IncKurt, "Tables/Table06_Include_Kurt_Lower_Panel.csv")

# Moments Only ----------------------------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=0) #Make.Reg.Data.Moments
myRegData.plus1 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=1)
myRegData.plus2 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=2)
myRegData.plus3 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + SD + Skew + Kurtosis + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData.clean <- subset(iRegData, select = -c(VIX))
  
  
  if (sum(!apply(iRegData.clean, 1, anyNA)) > 120) {
    print(iID)
    
    #iRegData.clean <- subset(iRegData.clean, Recession == 0) # Choose between 0 and 1
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(SD)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval", "SD", "SD_pval", "Skew", "Skew_pval","Kurtosis", "Kurtosis_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","SD_pval","Skew_pval","Kurtosis_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Moment <- matrix(nrow=9,ncol=8)
Res.ind.Summary.Moment[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Moment[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.Moment[3,] <- Cal.Summary(Res.ind[,"SD"], Res.ind[,"SD_pval"])
Res.ind.Summary.Moment[4,] <- Cal.Summary(Res.ind[,"Skew"], Res.ind[,"Skew_pval"])
Res.ind.Summary.Moment[5,] <- Cal.Summary(Res.ind[,"Kurtosis"], Res.ind[,"Kurtosis_pval"])
Res.ind.Summary.Moment[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.Moment[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Moment[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Moment[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Moment) <- c("Intercept","CPI", "SD","Skew","Kurtosis", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Moment) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Moment, "Tables/Table07_Moments_Only.csv")

# High/Low Economic Uncertainty -----------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=0) #Make.Reg.Data.Moments
myRegData.plus1 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=1)
myRegData.plus2 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=2)
myRegData.plus3 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

# High Economic Uncertainty 
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData <- subset(iRegData, select = -c(VIX))
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, Recession == 1) # Choose between 0 and 1
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.HighEU <- matrix(nrow=9,ncol=8)
Res.ind.Summary.HighEU[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.HighEU[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.HighEU[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.HighEU[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.HighEU[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.HighEU[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.HighEU[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.HighEU[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.HighEU[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.HighEU) <- c("Intercept","CPI", "SD","Skew","Kurtosis", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.HighEU) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.HighEU, "Tables/Table08_HighEU_Uppel_Panel.csv")

# Low Economic Uncertainty 
Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData <- subset(iRegData, select = -c(VIX))
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, Recession == 0) # Choose between 0 and 1
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.LowEU <- matrix(nrow=9,ncol=8)
Res.ind.Summary.LowEU[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.LowEU[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.LowEU[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.LowEU[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.LowEU[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.LowEU[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.LowEU[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.LowEU[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.LowEU[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.LowEU) <- c("Intercept","CPI", "SD","Skew","Kurtosis", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.LowEU) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.LowEU, "Tables/Table08_LowEU_Lower_Panel.csv")

# High/Low Economic Uncertainty for the Moments-Based Model -------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=0) #Make.Reg.Data.Moments
myRegData.plus1 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=1)
myRegData.plus2 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=2)
myRegData.plus3 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

# High Economic Uncertainty 
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + SD + Skew + Kurtosis + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData <- subset(iRegData, select = -c(VIX))
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, Recession == 1) # Choose between 0 and 1
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval", "SD", "SD_pval", "Skew", "Skew_pval","Kurtosis", "Kurtosis_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","SD_pval","Skew_pval","Kurtosis_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Moment.HighEU <- matrix(nrow=9,ncol=8)
Res.ind.Summary.Moment.HighEU[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Moment.HighEU[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.Moment.HighEU[3,] <- Cal.Summary(Res.ind[,"SD"], Res.ind[,"SD_pval"])
Res.ind.Summary.Moment.HighEU[4,] <- Cal.Summary(Res.ind[,"Skew"], Res.ind[,"Skew_pval"])
Res.ind.Summary.Moment.HighEU[5,] <- Cal.Summary(Res.ind[,"Kurtosis"], Res.ind[,"Kurtosis_pval"])
Res.ind.Summary.Moment.HighEU[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.Moment.HighEU[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Moment.HighEU[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Moment.HighEU[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Moment.HighEU) <- c("Intercept","CPI", "SD","Skew","Kurtosis", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Moment.HighEU) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Moment.HighEU, "Tables/Table09_Moments_HighEU_Upper_Panel.csv")

# Low Economic Uncertainty 
Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData <- subset(iRegData, select = -c(VIX))
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, Recession == 0) # Choose between 0 and 1
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval", "SD", "SD_pval", "Skew", "Skew_pval","Kurtosis", "Kurtosis_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","SD_pval","Skew_pval","Kurtosis_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Moment.LowEU <- matrix(nrow=9,ncol=8)
Res.ind.Summary.Moment.LowEU[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Moment.LowEU[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.Moment.LowEU[3,] <- Cal.Summary(Res.ind[,"SD"], Res.ind[,"SD_pval"])
Res.ind.Summary.Moment.LowEU[4,] <- Cal.Summary(Res.ind[,"Skew"], Res.ind[,"Skew_pval"])
Res.ind.Summary.Moment.LowEU[5,] <- Cal.Summary(Res.ind[,"Kurtosis"], Res.ind[,"Kurtosis_pval"])
Res.ind.Summary.Moment.LowEU[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.Moment.LowEU[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Moment.LowEU[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Moment.LowEU[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Moment.LowEU) <- c("Intercept","CPI", "SD","Skew","Kurtosis", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Moment.LowEU) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Moment.LowEU, "Tables/Table09_Moments_LowEU_Lower_Panel.csv")

# Encompassing J-test ----------------------------------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=0) #Make.Reg.Data.Moments
myRegData.plus1 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=1)
myRegData.plus2 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=2)
myRegData.plus3 <- Make.Reg.Data.Moments.External(CPI, UNEMP, Scores, Moments, External, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

# Full-Sample
myFormula.FPC <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0') #+ Mean
myFormula.Mom <- paste0('CPI.hplus0 ~ CPI.hplus1 + SD + Skew + Kurtosis + UNEMP.hplus0') #+ Mean

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData.clean <- subset(iRegData, select = -c(VIX))
  
  
  if (sum(!apply(iRegData.clean, 1, anyNA)) > 120) {
    print(iID)
    
    #iRegData.clean <- subset(iRegData, Recession == 1) # Choose between 0 and 1
    iRegData.clean <- subset(iRegData, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    
    iRes.FPC <- lm(myFormula.FPC, data=iRegData.clean)
    iRes.Mom <- lm(myFormula.Mom, data=iRegData.clean)
    iRes.jtest <- jtest(iRes.FPC, iRes.Mom)
    
    Res[[i]] <- c(iID, iRes.jtest$`Pr(>|t|)`[2], iRes.jtest$`Pr(>|t|)`[1])
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
Res.ind.Jtest.Full <- Count.Pairs(Res.ind[,2:3])
write.csv(Res.ind.Jtest.Full, "Tables/Table10_Jtest_Full_Upper_Panel.csv")

# High Economic Uncertainty
uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData.clean <- subset(iRegData, select = -c(VIX))
  
  
  if (sum(!apply(iRegData.clean, 1, anyNA)) > 120) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, Recession == 1) # Choose between 0 and 1
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    
    iRes.FPC <- lm(myFormula.FPC, data=iRegData.clean)
    iRes.Mom <- lm(myFormula.Mom, data=iRegData.clean)
    iRes.jtest <- jtest(iRes.FPC, iRes.Mom)
    
    Res[[i]] <- c(iID, iRes.jtest$`Pr(>|t|)`[2], iRes.jtest$`Pr(>|t|)`[1])
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
Res.ind.Jtest.HighEU <- Count.Pairs(Res.ind[,2:3])
write.csv(Res.ind.Jtest.HighEU, "Tables/Table10_Jtest_HighEU_Middle_Panel.csv")

# Low Economic Uncertainty
uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData.clean <- subset(iRegData, select = -c(VIX))
  
  
  if (sum(!apply(iRegData.clean, 1, anyNA)) > 120) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, Recession == 0) # Choose between 0 and 1
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    
    iRes.FPC <- lm(myFormula.FPC, data=iRegData.clean)
    iRes.Mom <- lm(myFormula.Mom, data=iRegData.clean)
    iRes.jtest <- jtest(iRes.FPC, iRes.Mom)
    
    Res[[i]] <- c(iID, iRes.jtest$`Pr(>|t|)`[2], iRes.jtest$`Pr(>|t|)`[1])
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
Res.ind.Jtest.LowEU <- Count.Pairs(Res.ind[,2:3])
write.csv(Res.ind.Jtest.LowEU, "Tables/Table10_Jtest_LowEU_Lower_Panel.csv")

# High/Low Inflation Times ----------------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=0)
myRegData.plus1 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=1)
myRegData.plus2 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=2)
myRegData.plus3 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

# High Inflation Times
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData <- subset(iRegData, select = -c(VIX))
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, Inflation > 2) # Choose between Inflation > 2 and Inflation <= 2; For inflation target, set it as 1 <= Inflation & Inflation < 3 or !(1 < Inflation & Inflation < 3)
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.HighInf <- matrix(nrow=9,ncol=8)
Res.ind.Summary.HighInf[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.HighInf[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.HighInf[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.HighInf[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.HighInf[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.HighInf[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.HighInf[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.HighInf[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.HighInf[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.HighInf) <- c("Intercept","CPI","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.HighInf) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.HighInf, "Tables/Table11_High_Inflation_Upper_Panel.csv")

# Lower Inflation Times
Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData <- subset(iRegData, select = -c(VIX))
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, Inflation <= 2) # Choose between Inflation > 2 and Inflation <= 2; For inflation target, set it as 1 <= Inflation & Inflation < 3 or !(1 < Inflation & Inflation < 3)
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.LowInf <- matrix(nrow=9,ncol=8)
Res.ind.Summary.LowInf[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.LowInf[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.LowInf[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.LowInf[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.LowInf[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.LowInf[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.LowInf[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.LowInf[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.LowInf[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.LowInf) <- c("Intercept","CPI","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.LowInf) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.LowInf, "Tables/Table11_Low_Inflation_Lower_Panel.csv")

# Inflation Target ------------------------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=0)
myRegData.plus1 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=1)
myRegData.plus2 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=2)
myRegData.plus3 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

# Close to Inflation Target
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData <- subset(iRegData, select = -c(VIX))
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, 1 <= Inflation & Inflation < 3) # Choose between Inflation > 2 and Inflation <= 2; For inflation target, set it as 1 <= Inflation & Inflation < 3 or !(1 < Inflation & Inflation < 3)
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Close_InfTar <- matrix(nrow=9,ncol=8)
Res.ind.Summary.Close_InfTar[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Close_InfTar[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.Close_InfTar[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.Close_InfTar[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.Close_InfTar[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.Close_InfTar[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.Close_InfTar[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Close_InfTar[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Close_InfTar[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Close_InfTar) <- c("Intercept","CPI","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Close_InfTar) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Close_InfTar, "Tables/Table12_Close_InfTar_Upper_Panel.csv")

# Lower Inflation Times
Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  iRegData <- subset(iRegData, select = -c(VIX))
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(1 < Inflation & Inflation < 3)) # Choose between Inflation > 2 and Inflation <= 2; For inflation target, set it as 1 <= Inflation & Inflation < 3 or !(1 < Inflation & Inflation < 3)
    iRegData.clean <- subset(iRegData.clean, select = -c(Recession))
    iRegData.clean <- subset(iRegData.clean, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Far_InfTar <- matrix(nrow=9,ncol=8)
Res.ind.Summary.Far_InfTar[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Far_InfTar[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.Far_InfTar[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.Far_InfTar[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.Far_InfTar[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.Far_InfTar[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.Far_InfTar[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Far_InfTar[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Far_InfTar[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Far_InfTar) <- c("Intercept","CPI","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Far_InfTar) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Far_InfTar, "Tables/Table12_Far_Inflation_Lower_Panel.csv")

# Inflation Gap ---------------------------------------------------------------------------------------------------------------------------------
CPI_Org <- CPI

# Make Inflation Gap
CPI$prev1 <- CPI$prev1 - CPI.long$CPI10
CPI$plus0 <- CPI$plus0 - CPI.long$CPI10
CPI$plus1 <- CPI$plus1 - CPI.long$CPI10
CPI$plus2 <- CPI$plus2 - CPI.long$CPI10
CPI$plus3 <- CPI$plus3 - CPI.long$CPI10
CPI$plus4 <- CPI$plus4 - CPI.long$CPI10

# Make Regression Data
myRegData.plus0 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=0)
myRegData.plus1 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=1)
myRegData.plus2 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=2)
myRegData.plus3 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <-  dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <-  dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Gap <- matrix(nrow=9,ncol=8)
Res.ind.Summary.Gap[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Gap[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.Gap[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.Gap[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.Gap[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.Gap[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.Gap[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Gap[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Gap[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Gap) <- c("Intercept","CPI","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Gap) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Gap, "Tables/Table13_Inflation_Gap.csv")

# Hybrid-PC -------------------------------------------------------------------------------------------------------------------------------------
CPI <- CPI_Org

myRegData.plus1 <- Make.Reg.Data.hybrid(CPI, UNEMP, Scores, h=1)
myRegData.plus2 <- Make.Reg.Data.hybrid(CPI, UNEMP, Scores, h=2)
myRegData.plus3 <- Make.Reg.Data.hybrid(CPI, UNEMP, Scores, h=3)
myRegData <- rbind(myRegData.plus1, myRegData.plus2, myRegData.plus3)

myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + CPI.hminus1 + score1 + score2 + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + CPI.hminus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 90) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(CPI.hminus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI_f","CPI_f_pval","CPI_b","CPI_b_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <-  dplyr::select(Res.DF,c("Intercept_pval","CPI_f_pval","CPI_b_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <-  dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Hybrid <- matrix(nrow=10,ncol=8)
Res.ind.Summary.Hybrid[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Hybrid[2,] <- Cal.Summary(Res.ind[,"CPI_f"], Res.ind[,"CPI_f_pval"])
Res.ind.Summary.Hybrid[3,] <- Cal.Summary(Res.ind[,"CPI_b"], Res.ind[,"CPI_b_pval"])
Res.ind.Summary.Hybrid[4,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.Hybrid[5,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.Hybrid[6,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.Hybrid[7,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.Hybrid[8,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Hybrid[9,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Hybrid[10,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Hybrid) <- c("Intercept","CPI_forward","CPI_backward","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Hybrid) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Hybrid, "Tables/Table14_Hybrid_PC.csv")

# Decompose Score2 ------------------------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.decompseS2(CPI, UNEMP, Scores, h=0)
myRegData.plus1 <- Make.Reg.Data.decompseS2(CPI, UNEMP, Scores, h=1)
myRegData.plus2 <- Make.Reg.Data.decompseS2(CPI, UNEMP, Scores, h=2)
myRegData.plus3 <- Make.Reg.Data.decompseS2(CPI, UNEMP, Scores, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2_pos + score2_neg  + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2_pos", "score2_pos_pval","score2_neg", "score2_neg_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pos_pval","score2_neg_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Decom2 <- matrix(nrow=10,ncol=8)
Res.ind.Summary.Decom2[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Decom2[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.Decom2[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.Decom2[4,] <- Cal.Summary(Res.ind[,"score2_pos"], Res.ind[,"score2_pos_pval"])
Res.ind.Summary.Decom2[5,] <- Cal.Summary(Res.ind[,"score2_neg"], Res.ind[,"score2_neg_pval"])
Res.ind.Summary.Decom2[6,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.Decom2[7,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP"])
Res.ind.Summary.Decom2[8,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Decom2[9,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Decom2[10,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Decom2) <- c("Intercept","CPI","Score_1","Score_2_Positive","Score_2_Negative", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Decom2) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Decom2, "Tables/Table15_Decompose_Score2.csv")

# Fixed Effects  --------------------------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=0)
myRegData.plus1 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=1)
myRegData.plus2 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=2)
myRegData.plus3 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=3)

myRegData.plus0$h <- 'h0' 
myRegData.plus1$h <- 'h1' 
myRegData.plus2$h <- 'h2' 
myRegData.plus3$h <- 'h3' 

myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)
myRegData$h <- as.factor(myRegData$h)

myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0 + h')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0 + h')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","h1","h1_pval","h2","h2_pval","h3","h3_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","h1_pval","h2_pval","h3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.FE <- matrix(nrow=12,ncol=8)
Res.ind.Summary.FE[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.FE[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.FE[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.FE[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.FE[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.FE[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.FE[7,] <- Cal.Summary(Res.ind[,"h1"], Res.ind[,"h1_pval"])
Res.ind.Summary.FE[8,] <- Cal.Summary(Res.ind[,"h2"], Res.ind[,"h2_pval"])
Res.ind.Summary.FE[9,] <- Cal.Summary(Res.ind[,"h3"], Res.ind[,"h3_pval"])
Res.ind.Summary.FE[10,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.FE[11,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.FE[12,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.FE) <- c("Intercept","CPI","Score_1","Score_2", "Score_3", "UNEMP", "FE_h1", "FE_h2", "FE_h3","F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.FE) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.FE, "Tables/Table16_Fixed_Effects.csv")

# Include VIX -----------------------------------------------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=0)
myRegData.plus1 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=1)
myRegData.plus2 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=2)
myRegData.plus3 <- Make.Reg.Data.External(CPI, UNEMP, Scores, External, h=3)
myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)

myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + VIX + score1 + score2 + score3 + UNEMP.hplus0') # CPI.hplus0 ~ CPI.hplus1 + VIX + score1 + score2 + score3 + UNEMP.hplus0
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1  + VIX + UNEMP.hplus0') # CPI.hplus0 ~ CPI.hplus1 + VIX + UNEMP.hplus0

## Individual Regression
uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) {
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(VIX)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval", "VIX", "VIX_pval", "score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","VIX_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.VIX <- matrix(nrow=10,ncol=8)
Res.ind.Summary.VIX[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.VIX[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.VIX[3,] <- Cal.Summary(Res.ind[,"VIX"], Res.ind[,"VIX_pval"])
Res.ind.Summary.VIX[4,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.VIX[5,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.VIX[6,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.VIX[7,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.VIX[8,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.VIX[9,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.VIX[10,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.VIX) <- c("Intercept","CPI", "VIX","Score_1","Score_2", "Score_3", "UNEMP","F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.VIX) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.VIX, "Tables/Table17_Include_VIX.csv")

# First Difference of Unemployment Rate Forecasts -----------------------------------------------------------------------------------------------
myRegData.plus0 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=0)
myRegData.plus1 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=1)
myRegData.plus2 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=2)
myRegData.plus3 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=3)

myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData$UNEMP.hplus0 <- c(NA, diff(iRegData$UNEMP.hplus0))
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    
    Res[[i]] <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    no.col <- length(Res[[i]])
  }
}

Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.FD <- matrix(nrow=9,ncol=8)
Res.ind.Summary.FD[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.FD[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.FD[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.FD[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.FD[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.FD[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.FD[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.FD[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.FD[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.FD) <- c("Intercept","CPI", "Score_1","Score_2", "Score_3", "FD_UNEMP","F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.FD) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.FD, "Tables/Table18_FD_Unemployment.csv")

# Remove Outliers at 5% -----------------------------------------------------------------------------------------------
Density <- list()
for (t in 1:nDate){
  t
  qtr <- UniDate[t]
  CPI.qtr <- subset(CPI, Date %in% qtr)
  
  for (h in 1:nH){
    h
    hname <- Hname[h]
    CPI.qtr.h <- CPI.qtr[[hname]]
    
    t.h.dens <- density(CPI.qtr.h, bw = bw.choice, from=min(CPI.qtr.h, na.rm = T), to=max(CPI.qtr.h, na.rm = T), n = npoint, na.rm=TRUE)
    t.h.from <- quantile(CPI.qtr.h, 0.05, na.rm = T) # 0.01
    t.h.to <- quantile(CPI.qtr.h, 0.95, na.rm = T) # 0.99
    # t.h.from <- -5
    # t.h.to <- 10
    t.h.dens <- density(CPI.qtr.h, bw = bw.choice, from=t.h.from, to=t.h.to, n = npoint, na.rm=TRUE)
    t.h.dens.dSup <- t.h.dens$x
    t.h.dens.norm <- t.h.dens %>% normalize.density() %>% RegulariseByAlpha(x=t.h.dens$x, alpha=alpha)
    
    Density[[hname]][['dens']][[t]] <- t.h.dens.norm
    Density[[hname]][['qd']][[t]]   <- dens2qd(dens=t.h.dens.norm, dSup = t.h.dens.dSup, qdSup =lqdSup)
    Density[[hname]][['lqd']][[t]]  <- dens2lqd(dens=t.h.dens.norm,dSup = t.h.dens.dSup, lqdSup=lqdSup)
    Density[[hname]][['dSup']][[t]] <- t.h.dens.dSup
    
    Density[[hname]][['Mean']][[t]] <- mean(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['Median']][[t]] <- median(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['SD']][[t]] <- sd(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['IQR']][[t]] <- IQR(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['Skew']][[t]] <- skewness(CPI.qtr.h, na.rm=TRUE)
    Density[[hname]][['Kurt']][[t]] <- kurtosis(CPI.qtr.h, na.rm=TRUE) + 3
    Density[[hname]][['Curv']][[t]] <- trapz(lqdSup, c(0,0,diff(diff(Density[[hname]][['lqd']][[t]])))^2)
  }
}

# Extract scores
lqd.plus1 <- matrix(unlist(Density[['plus1']]$lqd), nrow=nDate, ncol=npoint, byrow = T)
lqd.plus2 <- matrix(unlist(Density[['plus2']]$lqd), nrow=nDate, ncol=npoint, byrow = T)
lqd.plus3 <- matrix(unlist(Density[['plus3']]$lqd), nrow=nDate, ncol=npoint, byrow = T)
lqd.plus4 <- matrix(unlist(Density[['plus4']]$lqd), nrow=nDate, ncol=npoint, byrow = T)
lqd.pool  <- rbind(lqd.plus1, lqd.plus2, lqd.plus3, lqd.plus4)

LQD <- MakeFPCAInputs(tVec = lqdSup, yVec = lqd.pool)
FPCA.lqd <- fdapace::FPCA(Ly = LQD$Ly, Lt = LQD$Lt)

Mu <- FPCA.lqd$mu
PC <- FPCA.lqd$phi[,1:3]
scores.allh <- FPCA.lqd$xiEst[,1:3]

Scores <- list()
Scores$plus1 <- scores.allh[(nDate*0+1):(nDate*1),]
Scores$plus2 <- scores.allh[(nDate*1+1):(nDate*2),]
Scores$plus3 <- scores.allh[(nDate*2+1):(nDate*3),]
Scores$plus4 <- scores.allh[(nDate*3+1):(nDate*4),]

myRegData.plus0 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=0)
myRegData.plus1 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=1)
myRegData.plus2 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=2)
myRegData.plus3 <- Make.Reg.Data.h(CPI, UNEMP, Scores, h=3)

myRegData <- rbind(myRegData.plus0, myRegData.plus1, myRegData.plus2, myRegData.plus3)
myFormula <- paste0('CPI.hplus0 ~ CPI.hplus1 + score1 + score2 + score3 + UNEMP.hplus0')
myFormula.reduce <- paste0('CPI.hplus0 ~ CPI.hplus1 + UNEMP.hplus0')

# Individual Regression
uniID <- unique(myRegData$ID)
nID <- length(uniID)

Res <- list()
Res2 <- list()
for (i in 1:nID){
  i
  iID <- uniID[i]
  iRegData <- subset(myRegData, ID %in% iID)  
  
  
  if (sum(!apply(iRegData, 1, anyNA)) > 120) { # 30 forecasts per h
    print(iID)
    
    iRegData.clean <- subset(iRegData, !(is.na(CPI.hplus0)|is.na(CPI.hplus1)|is.na(score1)|is.na(UNEMP.hplus0)))
    iRes <- lm(myFormula, data=iRegData.clean) 
    iRes.summary <- summary(iRes)
    
    iRes.reduce <- lm(myFormula.reduce, data=iRegData.clean)
    iRes.Ftest <- anova(iRes.reduce, iRes)
    iRes.reduce.summary <- summary(iRes.reduce)
    
    Res[[i]]  <- c(iID, matrix(t(iRes.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.Ftest$`Pr(>F)`[2], iRes.summary$adj.r.squared, nobs(iRes))
    Res2[[i]] <- c(iID, matrix(t(iRes.reduce.summary$coefficients[,c('Estimate','Pr(>|t|)')]), nrow=1), iRes.reduce.summary$adj.r.squared, nobs(iRes.reduce))
    
    no.col <- length(Res[[i]])
    no.col2 <- length(Res2[[i]])
  }
}

# Results
Res.ind <- matrix(unlist(Res), ncol=no.col, byrow=T)
colnames(Res.ind) <- c("ID","Intercept","Intercept_pval","CPI","CPI_pval","score1", "score1_pval","score2", "score2_pval","score3", "score3_pval","UNEMP","UNEMP_pval","F_pval","Adj_R2", "No_Obs")
Res.DF <- as_tibble(Res.ind)
pval.DF <- dplyr::select(Res.DF,c("Intercept_pval","CPI_pval","score1_pval","score2_pval","score3_pval","UNEMP_pval","F_pval"))
diag.DF <- dplyr::select(Res.DF,c("Adj_R2","No_Obs"))

Res.ind.Summary.Outlier5 <- matrix(nrow=9,ncol=8)
Res.ind.Summary.Outlier5[1,] <- Cal.Summary(Res.ind[,"Intercept"], Res.ind[,"Intercept_pval"])
Res.ind.Summary.Outlier5[2,] <- Cal.Summary(Res.ind[,"CPI"], Res.ind[,"CPI_pval"])
Res.ind.Summary.Outlier5[3,] <- Cal.Summary(Res.ind[,"score1"], Res.ind[,"score1_pval"])
Res.ind.Summary.Outlier5[4,] <- Cal.Summary(Res.ind[,"score2"], Res.ind[,"score2_pval"])
Res.ind.Summary.Outlier5[5,] <- Cal.Summary(Res.ind[,"score3"], Res.ind[,"score3_pval"])
Res.ind.Summary.Outlier5[6,] <- Cal.Summary(Res.ind[,"UNEMP"], Res.ind[,"UNEMP_pval"])
Res.ind.Summary.Outlier5[7,6:8] <- Cal.Summary(Pval=Res.ind[,"F_pval"])
Res.ind.Summary.Outlier5[8,1:5] <- Cal.Summary(Res.ind[,"Adj_R2"])
Res.ind.Summary.Outlier5[9,1:5] <- Cal.Summary(Res.ind[,"No_Obs"])
row.names(Res.ind.Summary.Outlier5) <- c("Intercept","CPI","Score_1","Score_2", "Score_3", "UNEMP", "F-test", "Adj_R2", "No_Obs" )
colnames(Res.ind.Summary.Outlier5) <- c("mean", "sd","lq","median","uq","rej_1","rej_5","rej_10")
write.csv(Res.ind.Summary.Outlier5, "Tables/Table19_Outlier5.csv")