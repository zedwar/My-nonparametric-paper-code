library(fGarch) 
library(ghyp) 
##Data 
hh<-read.table("D:\dj_tick.csv",sep=",", skip=1, nrows=10000) 
pr<-hh[,4] plot(pr,type="l") 
ff<-(log(pr)) 
plot(ff,type="l") 
r<-diff(ff) 
plot(r,type="l")
Losses<-c(500) 
for(i in 500:9999) {
  if (r[i]<0) 
    Losses[i-499]=-r[i] 
  else 
    Losses[i-499]=0 
  }
tm<-hh[,3] 
ttt<-timeDate(paste(tm), format = "%H:%M:%S") 
dt<-diff(ttt) 
d<-as.numeric(dt)
rLosses<-c(500) 
eps=r/sqrt(d) 
for(i in 500:9999) { 
  if (r[i]<0) 
    rLosses[i-499]=-r[i] 
  else 
    rLosses[i-499]=0 
  } 
plot(rLosses[1:9499],type="h") 
##--------------------------------
## ACD parameters estimation
dur<-d[1:499]
lhood<-function(teta,dur) { 
  psi<-c(499,NA) 
  if (teta[1]/(1-teta[3])>0) { 
    psi[1]=teta[1]/(1-teta[3]) 
    for(i in 2:499) 
      psi[i]=teta[1]+teta[2]*dur[i-1]+teta[3]*psi[i-1] 
    lhood=sum(log(psi)+dur/psi) 
  } 
  else 
    lhood=Inf 
  }
predict<-function(dur,a,b,c){
  psi<-c(500,NA)
  psi[1]=a/(1-c) 
  for(i in 2:500)
    psi[i]=a+b*dur[i-1]+c*psi[i-1]
  predict=psi[500] 
  } 
##--------------------------------
## UHF-GARCH parameters estimation
ret<-r[1:499] 
eps<-ret/sqrt(dur)
Vlhood<-function(theta,eps,dur) { 
  sig<-c(499,NA) 
  if( theta[1]*(1-theta[2]-theta[3])> 0) { 
    sig[1]=sqrt(theta[1]/(1-theta[2]-theta[3]))
    for(i in 2:499) {
      if ((theta[1]+theta[2]*eps[i-1]^2+theta[3]*sig[i-1]^2+theta[4]/dur[i])>0) 
        sig[i]=sqrt(theta[1]+theta[2]*eps[i-1]^2+theta[3]*sig[i-1]^2+theta[4]/dur[i]) 
      else 
        sig[i]=Inf 
    } 
    Vlhood=sum(log(2*pi*sig^2)+(eps/sig)^2)/2 } 
  else 
    Vlhood=Inf }
predictVol<-function(vect,eps,d1) {
  sig<-c(500,NA) 
  sig[1]=sqrt(vect[1]/(1-vect[2]-vect[3])) 
  for(i in 2:500) 
    sig[i]=sqrt(vect[1]+vect[2]*eps[i-1]^2+vect[3]*sig[i-1]^2+vect[4]/d1[i]) 
  predictVol=sig[500] 
  }
Vol<-function(vect,eps,d1) {
  sig<-c(500,NA) 
  sig[1]=sqrt(vect[1]/(1-vect[2]-vect[3])) 
  for(i in 2:500) 
    sig[i]=sqrt(vect[1]+vect[2]*eps[i-1]^2+vect[3]*sig[i-1]^2+vect[4]/d1[i]) 
  Vol=sig 
  }
##----Shape of the distribution---------
shape<-function(f){
  x<-c(200000,NA) 
  y<-c(200000,NA) 
  x[1]=-0.0003 
  y[1]=dghyp(0.0001,fitgh) 
  for (i in 2:40) {
    k=-0.0001+i*0.00001 
    x[i]=k 
    y[i]=dghyp(k,fitgh) 
    } 
  plot(x,y,type="l") 
  }
shape(fitgh)
##My Quantile------------------------
f<-function(x) {
  f<-dghyp(x,fitgh)
  }
quant1<-function(q) { 
  k=0 
  p=1 
  a=-0.1 
  b=0 
  x=(a+b)/2 
  ind=0 
  while((abs(q-p)>1e-6)&&(k<=1000)) { 
    if (ind==0) {
      k1=k 
      k=k1+1 
      x=(a+b)/2 
      } 
    if(integrate(f, -0.1, x, stop.on.error = FALSE)$message=="OK") {
      ind=0 
      I<-integrate(f, -0.1, x, stop.on.error = FALSE) 
      p=I$value 
      if (q<p) { 
        b=x 
        side=-1 
        } 
      else { 
        a=x 
        side=1 
        } 
      } 
    else {
      x1=x 
      if (side==1) 
        x=x1+0.0001 
      else 
        x=x1-0.0001 
      ind=1 
      } 
    } 
  print(k) 
  print(p) 
  quant1<-x 
  } 
q<-quant1(0.01) 
q 
##----------------------------------
square<-function(f) { 
  I<-integrate(f, -0.1, 0.1) 
  square<-I$value 
  } 
print(square(fitgh)) 
##-----------------------
## Quantity of overlaps
Calc=function(var, n1, n2)
  { 
  num<-c(501,NA) 
  k=0 
  d=0 
  for (i in n1:n2) 
    if (var[i]<rLosses[i]) {
      k=d+1 
      d=k 
      num[k]=i 
      } 
  print(k) 
  print(num) 
  Calc=k 
  } 
##------------------------------------------
##------VaR----
tvar<-c(10000,NA) 
for (i in 499:9999) { 
  j=i-498 
  w<-d[j:i]
  r499<-r[j:i] 
  eps499<-r499/sqrt(w) 
  sg0=sum(eps499^2)/499
  teta0<-c(sg0,1/4,1/4) 
  ##
teta0<-c(1/2,1/4,1/4) 
opt<-nlminb(start = teta0, lhood, d = w) 
vect1<-c(opt$par[1],opt$par[2],opt$par[3])
d500<-w 
d500[500]=predict(w,vect1[1],vect1[2],vect1[3]) 
theta0<-c(sg0,3/16,3/16,2/16) 
##
theta0<-c(1/2,1/6,1/6,1/6) 
opt1<-nlminb(start = theta0, Vlhood, eps=eps499, d = w, lower=1e-10, upper=5) 
vect2<-c(opt1$par[1],opt1$par[2],opt1$par[3],opt1$par[4]) 
sig=Vol(vect2,eps499,d500) 
sigma=sig*sqrt(d500)
L=eps499/sigma[1:499] 
NIGfit<-fit.ghypuv(L) 
fitgh<-transform('NIGfit',0,sigma[500]) 
shape(fitgh) 
lambda<-quant1(0.01) 
tvar[j]=-lambda; 
} 
tvar[1:9501] 
plot(rLosses[1:9501],type="h") 
lines(tvar[1:9501], type="l", col="blue") 
Calc(tvar,1,9500) 
##---------------
