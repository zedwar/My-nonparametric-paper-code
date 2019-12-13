## Loading data points from the file-----------------------------------------------
library(MASS) 
library(fitdistrplus) 
library(ghyp)
z1<- read.xlsx("C:\\Users\\jimmy\\Desktop\\Data\\new.xlsx",1,stringsAsFactors=F) 
b1<- z1[ ,5] 
b1<- gsub(',','',b1)
b1 <- as.numeric(b1)
print(b1)
## Calculating log-returns ---------------------------------------------------------
x1<-diff(log(b1)) 
print(x1) 
plot(x1, type="h")
## Finding excessive loses sequence for Backtesting and obtaining it’s graphic ----
loses=c(length(x1), NA) 
for (i in 1:length(x1)) {
  if (x1[i]<0) 
    loses[i]=-x1[i] 
  else 
    loses[i]=0 
  } 
plot(loses, type="l")
## A function that predicts volatility one step forward ---------------------------
sigma=function(k,x,t) { 
  d=0 
  vec<-x[1:t] 
  vec1<-log(vec^2) 
  for(i in 1:k) 
    if(is.finite(vec1[t-i])) 
      d=d+vec1[t-i] 
  sigma=d/k }
## A function that calculates volatility for distribution’s parameters estimation -
sigmad=function(k,x,t) { 
  d=0 
  vec<-x[1:t] 
  vec1<-log(vec^2) 
  for(i in 0:k-1) 
    if(is.finite(vec1[t-i])) 
      d=d+vec1[t-i] 
  sigma=d/k 
  }
## Calculating volatility using window 36 to find VaR -----------------------------
volatility36<-c(913,NA) 
k=36 
{ for (i in 37:913) 
  volatility36[i-36]=sqrt(exp(sigma(k,x1,i))) 
  } 
volatility36
## Calculating volatility using window 36 to estimate distribution’s parameters ---
volatilityd36<-c(913,NA) 
k=36 
{ for (i in 36:913) 
  volatilityd36[i-35]=sqrt(exp(sigmad(k,x1,i))) 
  } 
volatilityd36
## Calculating VaR 1% using GENERALIZED HYPERBOLIC DISTRIBUTION -------------------
var36GHYP<-c(1877,NA) 
for (i in 534:913) { 
  j=i-498 
  a=j-35 
  b=i-35 
  c=i-34 
  window<-x1[j:i] 
  v1<-volatilityd36[a:b] 
  L=window/v1
  GHYPfit<-fit.ghypuv(L) 
  fitgh<-transform('GHYPfit',0,volatility36[c]) 
  lambda=qghyp(0.01,fitgh) 
  var36GHYP[i+1]=mean(fitgh)-lambda; 
  }
for (i in 535:913) { 
  if(is.na(var36GHYP[i])==TRUE) 
    var36GHYP[i]=var36GHYP[i-1] 
  }
## Calculating the number of overlaps and finding points where they occurred ------
Calc=function(varf) {
  num<-c(1342,NA) 
  k=0 
  d=0 
  for (i in 1:1342) 
    if (varf[i]<=Loses[i]) { 
      k=d+1 
      d=k 
      num[k]=i 
    } 
  print(k) 
  print(num) 
  Calc=k 
  }
## Drawing graphics of losses and VaR 1% ------------------------------------------
var36GHYPplot<-c(1342,NA) 
var36GHYPplot<-var36GHYP[535:913] 
Losses<-c(1342,NA) 
Losses<-loses[534:913] 
Losses[652]=0.1 
Loses<-loses[534:913] 
plot(Losses, type="h", col="blue") 
lines(Loses, type="h") 
lines(var36GHYPplot, type="l",col="red") 
Calc(var36GHYPplot)
## To calculate VaR 5% using GENERALIZED HYPERBOLIC DISTRIBUTION we follow the same steps that were shown above for calculating VaR 1% but the line "lambda=qghyp(0.01,fitgh)" should be changed to "lambda=qghyp(0.05,fitgh)"
## To calculate VaR using other distributions we take the same procedure as shown above but with some changes ----------------------------------------------------
## VaR 1% using HYPERBOLIC DISTRIBUTION -------------------------------------------... 
var36HYP<-c(1877,NA) 
for (i in 534:1877) {
  j=i-498 
  a=j-35 
  b=i-35 
  c=i-34 
  window<-x1[j:i] 
  v1<-volatilityd36[a:b] 
  L=window/v1 
  HYPfit<-fit.hypuv(L) 
  fith<-transform('HYPfit',0,volatility36[c]) 
  lambda=qghyp(0.01,fith) 
  var36HYP[i+1]=mean(fith)-lambda; 
  } 
## VaR 5% using HYPERBOLIC DISTRIBUTION -------------------------------------------
lambda=qghyp(0.01,fith) -----> lambda=qghyp(0.05,fith)
## VaR 1% using NORMAL INVERSE GAUSSIAN -------------------------------------------... 
var36NIG<-c(1877,NA) 
for (i in 534:1877) { 
  j=i-498 
  a=j-35 
  b=i-35 
  c=i-34 
  window<-x1[j:i] 
  v1<-volatilityd36[a:b] 
  L=window/v1 
  NIGfit<-fit.NIGuv(L) 
  fitnig<-transform('NIGfit',0,volatility36[c]) 
  lambda=qghyp(0.01,fitnig) 
  var36NIG[i+1]=mean(fitnig)-lambda; 
  } 
## VaR 5% using NORMAL INVERSE GAUSSIAN -------------------------------------------
lambda=qghyp(0.01,fitnig)-----> lambda=qghyp(0.05,fitnig)
## VaR 1% using NORMAL GAUSSIAN DISTRIBUTION --------------------------------------... 
var36GAUSS<-c(1877,NA) 
for (i in 534:1877) { 
  j=i-498 
  a=j-35 
  b=i-35 
  c=i-34 
  window<-x1[j:i] 
  v1<-volatilityd36[a:b] 
  L=window/v1 
  GAUSSfit<-fit.gaussuv(L) 
  par<-coef(GAUSSfit) 
  fitgauss<-gauss(par$mu, volatility36[c]*par$sigma)
  lambda=qghyp(0.01,fitgauss) var36GAUSS[i+1]=mean(fitgauss)-lambda; } ...
## VaR 5% using NORMAL GAUSSIAN DISTRIBUTION --------------------------------------
  lambda=qghyp(0.01,fitgauss)-----> lambda=qghyp(0.05,fitgauss)
## VaR 1% using Variance Gamma ----------------------------------------------------... 
ar36VG<-c(1877,NA) 
for (i in 534:1877) { 
  j=i-498 
  a=j-35 
  b=i-35 
  c=i-34 
  window<-x1[j:i] 
  v1<-volatilityd36[a:b] 
  L=window/v1 
  VGfit<-fit.VGuv(L) 
  fitvg<-transform('VGfit',0,volatility36[c]) 
  lambda=qghyp(0.01,fitvg) 
  var36VG[i+1]=mean(fitvg)-lambda; 
  }
## VaR 5% using Variance Gamma ----------------------------------------------------
lambda=qghyp(0.01,fitvg)-----> lambda=qghyp(0.05,fitvg)
## VaR 1% using STUDENT -----------------------------------------------------------... 
var36St<-c(1877,NA) 
for (i in 534:1877) { 
  j=i-498 
  a=j-35 
  b=i-35 
  c=i-34 
  window<-x1[j:i] 
  v1<-volatilityd36[a:b] 
  L=window/v1 
  Stfit<-fit.tuv(L) 
  fitst<-transform('Stfit',0,volatility36[c]) 
  lambda=qghyp(0.01,fitst) 
  var36St[i+1]=mean(fitst)-lambda;
  } 
## VaR 5% using STUDENT -----------------------------------------------------------
lambda=qghyp(0.01,fitst)-----> lambda=qghyp(0.05,fitst)

