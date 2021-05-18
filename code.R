#+eval=FALSE
#######################################################################################
# ECOM90022 Replication Task Mitchell Briggs 1048924                                  #
#                                                                                     #
# This code is a replication of Gabauer (2020), 'Volatility impulse response analysis #
# for DCC-GARCH models: The role of volatility transmission mechanisms'. Journal      #
# of Forecasting, 39(5), 788-796.                                                     # 
#                                                                                     #
# The original code can be found at https://github.com/GabauerDavid                   #
#######################################################################################

#Set working directory
setwd("C:/Users/Mitch/OneDrive/Desktop/UniMelb/2021/T1/ECOM90022 Research Methods/Replication Task")
options(digits=3)
options(scipen=999)
#Load required packages
library(CADFtest)
library(rmgarch)
library(rugarch)
library(tseries)
library(quantmod)
library(moments)
library(ggplot2)
library(stargazer)
library(parallel)
library(rmarkdown)

#Read in the data
dt <- read.csv("fx_gab.csv")

date <- dt$Date[2:4203]
as.Date(as.character(date))
#Generate time series objects - just for ts plots
et <- ts(dt$EUR, start=c(2002, 1), end=c(2018, 9), 262)
ft <- ts(dt$CHF, start=c(2002, 1), end=c(2018, 9), 262)
pt <- ts(dt$GBP, start=c(2002, 1), end=c(2018, 9), 262)
yt <- ts(dt$JPY, start=c(2002, 1), end=c(2018, 9), 262)

#Plot each series
par(mfrow=c(2, 2))
et1 <- ts.plot(et,
        main="EUR/USD Spot Rate",
        gpars = list(xlab="Year", ylab="EUR",
                     col=c("blue")))
ft1 <- ts.plot(ft,
        main="CHF/USD Spot Rate",
        gpars = list(xlab="Year",ylab="CHF",
                     col=c("blue")))
pt1 <- ts.plot(pt,
        main="GBP/USD Spot Rate",
        gpars = list(xlab="Year", ylab="GBP",
                     col=c("blue")))
yt1 <- ts.plot(yt,
        main="JPY/USD Spot Rate",
        gpars = list(xlab="Year", ylab="JPY",
                     col=c("blue")))

#Name each FX rate series
eur <- dt$EUR
franc <- dt$CHF
pound <- dt$GBP
yen <- dt$JPY

#Check for unit root in levels for each FX rate series
df_eur <-CADFtest(eur, type="drift", criterion="MAIC", max.lag.y=8)
etdf <- df_eur$p.value
df_franc <-CADFtest(franc, type="drift", criterion="MAIC", max.lag.y=8)
ftdf <- df_franc$p.value
df_pound <-CADFtest(pound, type="drift", criterion="MAIC", max.lag.y=8)
ptdf <- df_pound$p.value
df_yen <-CADFtest(yen, type="drift", criterion="MAIC", max.lag.y=8)
ytdf <- df_yen$p.value
df <- cbind(etdf, ftdf, ptdf, ytdf)
stargazer(df,
          digits=3,
          title="DF tests",
          colnames = TRUE,
          rownames = TRUE,
          type = "html",
          out="dftest.html",
          out.header = FALSE)

#Jarque & Bera (JB) normality test in levels
jb_eur <- jarque.test(eur)
jb_franc <- jarque.test(franc)
jb_pound <- jarque.test(pound)        
jb_yen <- jarque.test(yen) 
jb <- cbind(jb_eur$statistic, jb_franc$statistic, jb_pound$statistic, jb_yen$statistic)
stargazer(jb,
          digits=3,
          title="Jarque-Berra Normality Test",
          type = "html",
          colnames = TRUE,
          out="jbtest.html",
          out.header = FALSE)

#Compute returns for each FX rate series
l_eur <- diff(log(eur))*100
l_franc <- diff(log(franc))*100
l_pound <- diff(log(pound))*100
l_yen <- diff(log(yen))*100
#Unit root tests for returns data
df_eur2 <-CADFtest(l_eur, type="drift", criterion="MAIC", max.lag.y=8)
df_eur2$p.value
df_franc2 <-CADFtest(l_franc, type="drift", criterion="MAIC", max.lag.y=8)
df_franc2$p.value
df_pound2 <-CADFtest(l_pound, type="drift", criterion="MAIC", max.lag.y=8)
df_pound2$p.value
df_yen2 <-CADFtest(l_yen, type="drift", criterion="MAIC", max.lag.y=8)
df_yen2$p.value

l_et <- ts(l_eur, start=c(2002, 1), end=c(2018, 9), 262)
l_ft <- ts(l_franc, start=c(2002, 1), end=c(2018, 9), 262)
l_pt <- ts(l_pound, start=c(2002, 1), end=c(2018, 9), 262)
l_yt<- ts(l_yen, start=c(2002, 1), end=c(2018, 9), 262)

#Plot daily FX returns
par(mfrow=c(2, 2))
ts.plot(l_et,
        main="EUR/USD Daily Returns",
        gpars = list(xlab="Year", ylab="EUR",
                     col=c("blue")))
ts.plot(l_ft,
        main="CHF/USD Daily Returns",
        gpars = list(xlab="Year",ylab="CHF",
                     col=c("blue")))
ts.plot(l_pt,
        main="GBP/USD Daily Returns",
        gpars = list(xlab="Year", ylab="GBP",
                     col=c("blue")))
ts.plot(l_yt,
        main="JPY/USD Daily Returns",
        gpars = list(xlab="Year", ylab="JPY",
                     col=c("blue")))

#Combine log returns data in a single frame
log_data <- cbind(l_eur, l_franc, l_pound, l_yen)
as.data.frame(log_data)
log_data

#Compute unconditional correlations between FX rates
cor_fx <- cor(log_data)
colnames(cor_fx) <- c("EUR", "CHF", "GBP", "JPY")
rownames(cor_fx) <- c("EUR", "CHF", "GBP", "JPY")
stargazer(cor_fx,
          digits=3,
          title="Unconditional Correlation",
          type = "html",
          colnames = TRUE,
          out="corr.html",
          out.header = FALSE)

#Construct table of means for each FX rate
table_1 <- rbind(apply(log_data,2,mean),
             apply(log_data,2,var),
             apply(log_data,2,skewness),
             apply(log_data,2,kurtosis))
rownames(table_1) <- c("Mean","Variance","Skewness","Kurtosis")
colnames(table_1) <- c("EUR", "CHF", "GBP", "JPY")
round(table_1,3)
stargazer(table_1,
          digits=3,
          title="Decriptive Statistics",
          type = "html",
          colnames = TRUE,
          out="desstat.html",
          out.header = FALSE)


##Assume same GARCH(1,1) spec for each FX series
Y = log_data
k = ncol(Y)
NAMES = colnames(Y)
DATE = as.Date(as.character(dt[,1]))
date = DATE[-1]
t = nrow(Y)
#Specify univariate GARCH
garch11 = ugarchspec(mean.model=list(armaOrder=c(0,0)), 
                          variance.model=list(garchOrder=c(1,1), model="sGARCH"),
                          distribution.model="norm")
dcc.garch11.spec = dccspec(uspec=multispec(replicate(4, garch11)), 
                           dccOrder=c(1,1), distribution="mvnorm")
dcc.fit = dccfit(dcc.garch11.spec, data=Y)

#DCC test
dcc_test <- DCCtest(Y, garchOrder=c(1,1), solver="solnp", n.lags=8)
dcc_test

#Save estimation results
dcc_est <- dcc.fit@mfit$matcoef
stargazer(dcc_est,
          digits=3,
          title="DCC-GARCH(1,1) Estimation",
          type = "html",
          colnames = TRUE,
          out.header = FALSE)

#Plot dynamics of DCC-GARCH 
cor1 <- rcor(dcc.fit)  
voc1 <- rcov(dcc.fit)
cor1
par(mfrow=c(3,1))
ts.plot(cor1[1,2,], col=c("blue"))
ts.plot(cor1[1,3,], col=c("blue"))
ts.plot(cor1[1,4,], col=c("blue"))

#GFVED
DCA = function(CV, digit=2){
  k = dim(CV)[1]
  CT = apply(CV,1:2,mean)*100 
  OWN = diag(diag(CT))
  TO = colSums(CT-OWN)
  FROM = rowSums(CT-OWN)
  NET = TO-FROM
  TCI = mean(TO)
  NPSO = CT-t(CT)
  NPDC = rowSums(NPSO>0)
  table = format(round(cbind(CT,FROM),digit),nsmall=digit)
  to = c(format(round(c(TO,sum(TO)),digit),nsmall=digit))
  net = c(format(round(c(NET, TCI),digit),nsmall=digit))
  inc = c(format(round(colSums(CT), digit),nsmall=digit), "TCI")
  npdc = c(format(round(NPDC,digit),nsmall=digit), TCI*k/(k-1))
  TABLE = rbind(table,to,inc,net,npdc)
  colnames(TABLE) = c(rownames(CV),"FROM others")
  rownames(TABLE) = c(rownames(CV),"TO others","Inc. own","NET","NPDC")
  PCI = matrix(NA, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      PCI[i,j] = 2*(CT[i,j]+CT[j,i])/(CT[i,i]+CT[i,j]+CT[j,i]+CT[j,j])
    }
  }
  return = list(CT=CT,TCI=TCI,TCI_corrected=TCI*k/(k-1),PCI=PCI,TO=TO,FROM=FROM,NET=NET,NPSO=NPSO,NPDC=NPDC,TABLE=TABLE)
}
VGFEVD = function(dcc.fit, h=100, standardize=FALSE) {
  NAMES = dcc.fit@model$modeldata$asset.names
  R = rcor(dcc.fit)
  R.bar = apply(R,1:2,mean)
  Q.bar = dcc.fit@mfit$Qbar
  H = rcov(dcc.fit)
  t = dim(H)[3]
  #alpha
  alpha = array(0,c(k,k,h))
  alpha[,,1] = diag(dcc.fit@mfit$matcoef[c(seq(3,(4*k),4)),1])
  #beta
  beta = diag(dcc.fit@mfit$matcoef[c(seq(4,(4*k),4)),1])
  #ALPHA
  ALPHA = dcc.fit@mfit$matcoef[(4*k+1),1]
  #BETA
  BETA = dcc.fit@mfit$matcoef[(4*k+2),1]
  
  H.hat = array(0,c(k,k,h+1))
  GVIRF = H.hat.shock = H.hat.no_shock = array(0,c(k,k,t,h+1))
  e = diag(k)
  for (i in 1:t) {
    H.hat[,,1] = H[,,i]
    Q.hat = H.hat
    Q.hat[,,1] = dcc.fit@mfit$Q[[i]]
    for (j in 1:h) {
      H.hat[,,j+1] = (alpha[,,j])%*%e^2 + beta%*%H.hat[,,j]
      D = diag(diag(H.hat[,,j+1])^0.5)
      u = D%*%e
      if (j==1) {
        Q.hat[,,2] = (1-ALPHA-BETA)*Q.bar + ALPHA*crossprod(u) + BETA*H.hat[,,1]
      } else {
        Q.hat[,,j+1] = (1-ALPHA-BETA)*Q.bar + (ALPHA+BETA)*Q.hat[,,j]
      }
      R.hat = diag(1/(diag(Q.hat[,,j+1])^0.5))%*%Q.hat[,,j+1]%*%(diag(1/diag(Q.hat[,,j+1])^0.5))
      H.hat[,,j+1] = D%*%R.hat%*%D
    }
    H.hat.shock[,,i,] = H.hat
  }
  if (standardize) {
    e =0*diag(k)
    for (i in 1:t) {
      H.hat[,,1] = H[,,i]
      Q.hat = H.hat
      Q.hat[,,1] = dcc.fit@mfit$Q[[i]]
      for (j in 1:h) {
        H.hat[,,j+1] = beta%*%H.hat[,,j]
        D = diag(diag(H.hat[,,j+1])^0.5)
        if (j==1) {
          Q.hat[,,2] = (1-ALPHA-BETA)*Q.bar + BETA*H.hat[,,1]
        } else {
          Q.hat[,,j+1] = (1-ALPHA-BETA)*Q.bar+(ALPHA+BETA)*Q.hat[,,j]
        }
        R.hat = diag(1/(diag(Q.hat[,,j+1])^0.5))%*%Q.hat[,,j+1]%*%(diag(1/diag(Q.hat[,,j+1])^0.5))
        H.hat[,,j+1] = D%*%R.hat%*%D
      }
      H.hat.no_shock[,,i,] = H.hat
    }
    #for (i in 1:t) {
    # GVIRF[,,i,] = H.hat.shock[,,i,] - H.hat.no_shock[,,i,]
    #}
  }# else {
  for (i in 1:t) {
    GVIRF[,,i,] = H.hat.shock[,,i,] - H.hat.no_shock[,,i,]
  }
  #GVIRF = H.hat.shock
  #}
  GVFEVD = array(NA, c(k,k,t))
  rownames(GVFEVD) = colnames(GVFEVD)=NAMES
  for (i in 1:t) {
    num = apply(GVIRF[,,i,]^2,1:2,sum)
    den = c(apply(num,1,sum))
    fevd = t(num)/den
    GVFEVD[,,i] = (fevd/apply(fevd, 1, sum))
  }
  return = list(GVFEVD=GVFEVD, GVIRF=GVIRF)
}

nfore <- 100
dcc_gfevd <- VGFEVD(dcc.fit,nfore)
virf <- dcc_gfevd$GVIRF

#Connectedness Analysis
#Static connectedness table
table_2 <- (DCA(dcc_gfevd$GVFEVD)$TABLE)
table_2
stargazer(table_2,
          digits=3,
          title="Static Connectedness",
          type = "html",
          colnames = TRUE,
          out="static.html",
          out.header = FALSE)

net = matrix(NA,nrow=t,ncol=k)
to = matrix(NA,nrow=t,ncol=k)
from = matrix(NA,nrow=t,ncol=k)
total = matrix(NA,ncol=1,nrow=t)
ct = array(NA,c(k,k,t))
npso = array(NA,c(k,k,t))
for (i in 1:t) {
  dca = DCA(dcc_gfevd$GVFEVD[,,i])
  ct[,,i] = dca$CT
  total[i,] = dca$TCI_corrected
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
}


#Total volatility spillovers plot
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(total, type="l",xaxs="i",col="blue", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)


#Total directional 'TO' volatility spillovers plot
split <- 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(to[,i], xlab="",ylab="",type="l",xaxs="i",col="blue", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(0,ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
}

#Total directional 'FROM' volatility spillovers plot
par(oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(from[,i], xlab="",ylab="",type="l",xaxs="i",col="blue", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(0,ceiling(max(from))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
}

#Net 'TOTAL' directional spillovers plot
split <- 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(net[,i], xlab="",ylab="",type="l",xaxs="i",col="blue", las=1, main=paste("NET",colnames(Y)[i]),ylim=c(floor(min(net)), ceiling(max(net))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
}

#Net 'PAIRWISE' directional spillovers plot
kk <- k*(k-1)/2
lk <- ceiling(sqrt(2))
par(mfcol=c(ceiling(kk/lk),lk), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
  for (j in 1:k) {
    if (i<j) {
      plot(npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="blue", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
      grid(NA,NULL)
    }
  }
}

#End

