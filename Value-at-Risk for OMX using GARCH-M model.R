## Take 500 work data points 
z <- read.table("H:\omx1.txt", sep="", skip=4) 
y <- z[ ,5] 
b=c(length(y),NA) 
for (i in 1:length(y)) 
  b[i]=y[length(y)-i+1] 
x<-diff(log(b)) 
print(x)
## Backtesting 459 data points 
z <- read.table("H:\omx1.txt", sep= "", skip=4, nrows=500) 
y <- z[ ,2] 
x<-diff(log(y)) 
print(x)
z1 <- read.table("H:\omx1.txt", sep= "", skip=4) 
y1 <- z1[ ,2] 
x1<-diff(log(y1)) 
x2<-diff(y1)
## Excessive loses sequence for the Backtesting data 
loses=c(1378, NA) 
for (i in 500:1877) {
  if (x1[i]<0) 
    loses[i-499]=-x1[i] 
  else 
    loses[i-499]=0 
  } 
plot(loses, type="h")
## Excessive PRICE loses sequence for the Backtesting data 
Ploses=c(1378, NA) 
for (i in 500:1877) {
  if (x2[i]<0) 
    Ploses[i-499]=-x2[i] 
  else 
    Ploses[i-499]=0 
  } 
plot(Ploses, type="h") 
##-----------------------------
## Volatility GARCH sequence 
fitOmx = garchFit(~garch(1, 1), x) 
Interactive Plot: 
  plot(fitOmx, which=2)
##devolatilized returns 
OmxL=x/fit@sigma.t 
print(OmxL)
##GH fit 
fitted<-fit.ghypuv(OmxL, silent=TRUE) 
hist(Omxfit, ghyp.col="blue") 
qqghyp(Omxfit, gaussian=FALSE) 
##--------------------
## Quantity of overlaps
Calc=function(varf){ 
  num<-c(2008,NA) 
  k=0 
  d=0 
  for (i in 500:1877)
    if (varf[i]<=loses[i-499]) {
      k=d+1 
      d=k 
      num[k]=i 
      } 
  print(k) 
  print(num)
  Calc=k 
  } 
##----------------
## myGarch 1 step
mygarch=function(sigma1,L1,m,c,a,b) {
  mygarch=sqrt(c+a*sigma1*sigma1*(L1-m)*(L1-m)+ b*sigma1*sigma1) 
  } 
##---------------
## VaR 1:499
var<-c(1878,NA) 
for (i in 1:499) { 
  fitgh<-transform('fitted',0,fitOmx@sigma.t[i]) 
  lambda=qghyp(0.01,fitgh) 
  var[i]=-lambda; 
  } 
##----------
##(GH) Value-at-Risk
var<-c(1878,NA) 
Curve<-c(1378,NA) 
for (i in 499:1877) { 
  j=i-498 
  window<-x1[j:i] 
  g=garchFit(~garch(1, 1), window) 
  L=window/g@sigma.t 
  Ggfit<-fit.ghypuv(L) 
  sd=mygarch(g@sigma.t[499],L[499],g@fit$coef[1], g@fit$coef[2],g@fit$coef[3],g@fit$coef[4]) 
  Curve[j]=sd 
  fitgh<-transform('Ggfit',0,sd) 
  lambda=qghyp(0.01,fitgh) 
  var[i+1]=-lambda; 
  } 
plot(Curve) 
Calc(var) 
Rvar<-var[500:1877] 
Losses<-c(1879,NA) 
Losses<-loses[1:1378] 
Losses[1379]=0.15 
plot(Losses, type="h", col="white") 
lines(loses,type="h") 
##-------------------------
##Take Image out in the file
jpeg(filename="H:\QQplot for GH.jpg") 
qqghyp(fitted,gaussian=FALSE) 
dev.off() 
##----------
