## Loading data points from the file-----------------------------------------------
z1<- read.table("D:\dax.csv",sep=",", skip=2) 
y1<- z1[ ,5] 
b1=c(length(y1),NA) 
for (i in 1:length(y1)) 
  b1[i]=y1[length(y1)-i+1] 
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
plot(loses, type="h") 
lines(loses, type="h", col="blue")
## A function that predicts volatility one step forward ---------------------------
sigma=function(k,x,t) { 
  d=0 
  vec<-x[1:t] 
  vec1<-log(vec^2) 
  for(i in 1:k) 
    if(is.finite(vec1[t-i])) 
      d=d+vec1[t-i] 
  sigma=d/k 
  }
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
## Calculating volatility using window 35 to find VaR -----------------------------
volatility35<-c(1973,NA) 
k=35 
{ for (i in 36:2008) 
  volatility35[i-35]=sqrt(exp(sigma(k,x1,i))) 
  } 
volatility35
## Calculating volatility using window 35 to estimate distribution’s parameters ---
volatilityd35<-c(1973,NA) 
k=35 
{for (i in 35:2008) 
  volatilityd35[i-34]=sqrt(exp(sigmad(k,x1,i))) 
  } 
volatilityd35
## Calculating VaR 1% using GENERALIZED HYPERBOLIC DISTRIBUTION -------------------
var35GHYP<-c(2008,NA) 
for (i in 533:2008) { 
  j=i-498 
  a=j-34 
  b=i-34 
  c=i-33 
  window<-x1[j:i] 
  v1<-volatilityd35[a:b] 
  L=window/v1 
  GHYPfit<-fit.ghypuv(L) 
  fitgh<-transform('GHYPfit',0,volatility35[c]) 
  lambda=qghyp(0.01,fitgh) 
  var35GHYP[i+1]=mean(fitgh)-lambda; 
  }
## Calculating the number of overlaps and finding points where they occurred ------
Calc=function(varf) { 
  num<-c(2008,NA)
  k=0 
  d=0 
  for (i in 1:1474) 
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
var35GHYPplot<-c(1475,NA) 
var35GHYPplot<-var35GHYP[534:2007] 
Losses<-c(1475,NA) 
Losses<-loses[534:2007] 
Losses[1475]=0.18 
Loses<-loses[534:2007] 
plot(Losses, type="h", col="white") 
lines(Loses, type="h") 
lines(var35GHYPplot, type="l",col="red") 
Calc(var35GHYPplot) 
var35GHYPplot
## To calculate VaR 5% using GENERALIZED HYPERBOLIC DISTRIBUTION we follow the same steps that were shown above for calculating VaR 1% but the line "lambda=qghyp(0.01,fitgh)" should be changed to "lambda=qghyp(0.05,fitgh)"
## To calculate VaR using other distributions we take the same procedure as shown above but with following changes -----------------------------------------------
## VaR 1% using HYPERBOLIC DISTRIBUTION -------------------------------------------...
var35HYP<-c(2008,NA) 
for (i in 533:2008) { 
  j=i-498 
  a=j-34 
  b=i-34 
  c=i-33 
  window<-x1[j:i] 
  v1<-volatilityd35[a:b] 
  L=window/v1 
  HYPfit<-fit.hypuv(L) 
  fith<-transform('HYPfit',0,volatility35[c]) 
  lambda=qghyp(0.01,fith) 
  var35HYP[i+1]=mean(fith)-lambda; 
  }
## VaR 5% using HYPERBOLIC DISTRIBUTION -------------------------------------------
lambda=qghyp(0.01,fith) -----> lambda=qghyp(0.05,fith)
## VaR 1% using NORMAL INVERSE GAUSSIAN -------------------------------------------... 
var35NIG<-c(2008,NA) 
for (i in 533:2008) {
  j=i-498 
  a=j-34 
  b=i-34 
  c=i-33 
  window<-x1[j:i] 
  v1<-volatilityd35[a:b] 
  L=window/v1 NIGfit<-fit.NIGuv(L) 
  fitnig<-transform(‘NIGfit‘,0,volatility35[c]) 
  lambda=qghyp(0.01,fitnig) 
  var35NIG[i+1]=mean(fitnig)-lambda; 
  } 
## VaR 5% using NORMAL INVERSE GAUSSIAN -------------------------------------------
lambda=qghyp(0.01,finig) -----> lambda=qghyp(0.05,fitnig)
## VaR 1% using NORMAL GAUSSIAN DISTRIBUTION --------------------------------------... 
var35GAUSS<-c(2008,NA) 
for (i in 533:2008) { 
  j=i-498 
  a=j-34 
  b=i-34 
  c=i-33 
  window<-x1[j:i] 
  v1<-volatilityd35[a:b] 
  L=window/v1 
  par<-coef(GAUSSfit) 
  fitgauss<-gauss(par$mu, volatility35[c]*par$sigma) 
  lambda=qghyp(0.01,fitgauss) 
  var35GAUSS[i+1]=mean(fitgauss)-lambda; 
} 
## VaR 5% using NORMAL GAUSSIAN DISTRIBUTION --------------------------------------
lambda=qghyp(0.01,fitgauss)-----> lambda=qghyp(0.05,figauss)
## VaR 1% using Variance Gamma ----------------------------------------------------... 
var35VG<-c(2008,NA)
for (i in 533:2008) { 
  j=i-498 
  a=j-34 
  b=i-34 
  c=i-33 
  window<-x1[j:i] 
  v1<-volatilityd35[a:b] 
  L=window/v1 
  VGfit<-fit.VGuv(L) 
  fitvg<-transform(‘VGfit‘,0,volatility35[c]) 
  lambda=qghyp(0.01,fitvg) 
  var35VG[i+1]=mean(fitvg)-lambda; 
  } 
## VaR 5% using Variance Gamma ----------------------------------------------------
lambda=qghyp(0.01,fitvg)-----> lambda=qghyp(0.05,fitvg)
## VaR 1% using STUDENT -----------------------------------------------------------... 
var35St<-c(2008,NA) 
for (i in 533:2008) {
  j=i-498 
  a=j-34 
  b=i-34 
  c=i-33 
  window<-x1[j:i] 
  v1<-volatilityd35[a:b] 
  L=window/v1 
  Stfit<-fit.tuv(L) 
  fitst<-transform(‘Stfit‘,0,volatility35[c]) 
  lambda=qghyp(0.01,fitst) 
  var35St[i+1]=mean(fitst)-lambda; 
  } 
## VaR 5% using STUDENT -----------------------------------------------------------
lambda=qghyp(0.01,fitst) -----> lambda=qghyp(0.05,fitst)


