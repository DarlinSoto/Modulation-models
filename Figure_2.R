##########################
# AUTHOR: DARLIN SOTO    #
# DATE: OCTOBER 15, 2021 #
# EMAIL: dmsoto1@uc.cl   #
##########################

#This code reproduces Figure 2 in the paper Periodic Variable Stars Modulated by Time-Varying Parameters.

library(splines)

n <- 500
x <- sort(round(runif(n,0,1),5))

##############################################################################################################
##############################################################################################################
#SIMULATION 1: sinosoidal trend and amplitudes

aux1=2; aux2=9; aux3=6; aux4=4; aux5=7; w1=40*pi; w2=100*pi; lam_selec=c(50,1,2,10,1)
K=2
sigma2=2

mm=sin(aux1*pi*x)
g11=cos(aux2*pi*x)
g21=sin(aux3*pi*x)
g12=cos(aux4*pi*x)
g22=sin(aux5*pi*x)
C1=diag(cos(w1*x))
S1=diag(sin(w1*x))
C2=diag(cos(w2*x))
S2=diag(sin(w2*x))

mu=mm+C1%*%g11+S1%*%g21+C2%*%g12+S2%*%g22
plot(x,mu)

w_hat=c(w1,w2)
pordd=2
num_knots=30

set.seed(123)

coef_hat=matrix(0,ncol=SIM,nrow=((2*K+1)*(num_knots+3)))
y_hat=mu_hat=matrix(0,ncol=SIM,nrow=n)
sd_y_hat=matrix(0,ncol=SIM,nrow=n)
m_hat=sd_m_hat=matrix(0,ncol=SIM,nrow=n)
g11_hat=g21_hat=g12_hat=g22_hat=matrix(0,ncol=SIM,nrow=n)
sd_g11_hat=sd_g21_hat=sd_g12_hat=sd_g22_hat=matrix(0,ncol=SIM,nrow=n)
resid=matrix(0,ncol=SIM,nrow=n)
sigma2_hat=sigma2_hat2=rep(0,SIM)

SIM=200

for(m in 1:SIM){
  print(c(m,SIM))
  z1=rnorm(n,mean=0,sd=sqrt(sigma2))
  Y_t=mu+z1
  fit=psfit(x=x,xl=min(x)-0.0000001,xr=max(x)+0.0000001,y=Y_t,w1=w_hat,K=2,pord=pordd,
             ndx=num_knots,bdeg=3,lam=lam_selec)

  coef_hat[,m]=fit$theta
  y_hat[,m]=fit$f.hat
  mu_hat[,m]=fit$f.hat
  m_hat[,m]=fit$hat_m_g[,1]
  g11_hat[,m]=fit$hat_m_g[,2]
  g21_hat[,m]=fit$hat_m_g[,3]
  g12_hat[,m]=fit$hat_m_g[,4]
  g22_hat[,m]=fit$hat_m_g[,5]
  
  sigma2_hat[m]=fit$sigma2hat
  sigma2_hat2[m]=fit$sigma2hat2
  
  sd_y_hat[,m]=fit$sd_y
  sd_m_hat[,m]=fit$sd_m_g[,1]
  sd_g11_hat[,m]=fit$sd_m_g[,2]
  sd_g21_hat[,m]=fit$sd_m_g[,3]
  sd_g12_hat[,m]=fit$sd_m_g[,4]
  sd_g22_hat[,m]=fit$sd_m_g[,5]
}

#Response variable
y_hat_CI=apply(y_hat,1, quantile, probs=c(0.025,0.975))
y_hat_mean=apply(y_hat,1, mean)
sd_y_hat_mean=apply(sd_y_hat,1,mean)

mu_hat_CI=apply(mu_hat,1, quantile, probs=c(0.025,0.975))
mu_hat_mean=apply(mu_hat,1, mean)
sd_mu_hat_mean=sd_y_hat_mean

#trend and amplitudes
m_hat_CI=apply(m_hat,1, quantile, probs=c(0.025,0.975))
m_hat_mean=apply(m_hat,1, mean)
sd_m_hat_mean=apply(sd_m_hat,1,mean)
g11_hat_CI=apply(g11_hat,1, quantile, probs=c(0.025,0.975))
g11_hat_mean=apply(g11_hat,1, mean)
sd_g11_hat_mean=apply(sd_g11_hat,1, mean)
g21_hat_CI=apply(g21_hat,1, quantile, probs=c(0.025,0.975))
g21_hat_mean=apply(g21_hat,1, mean)
sd_g21_hat_mean=apply(sd_g21_hat,1, mean)
g12_hat_CI=apply(g12_hat,1, quantile, probs=c(0.025,0.975))
g12_hat_mean=apply(g12_hat,1, mean)
sd_g12_hat_mean=apply(sd_g12_hat,1, mean)
g22_hat_CI=apply(g22_hat,1, quantile, probs=c(0.025,0.975))
g22_hat_mean=apply(g22_hat,1, mean)
sd_g22_hat_mean=apply(sd_g22_hat,1, mean)

#Labels Y
Y_AR0_K2_SIN=Y_t
mu_AR0_K2_SIN_true=mu
mu_AR0_K2_SIN_esti=mu_hat_mean
mu_AR0_K2_SIN_sd=sd_mu_hat_mean
mu_AR0_K2_SIN_CI_Inf_1=mu_hat_CI[1,]
mu_AR0_K2_SIN_CI_Sup_1=mu_hat_CI[2,]

#Labels trend and amplitudes
m_AR0_K2_SIN_true=mm
m_AR0_K2_SIN_esti=m_hat_mean
m_AR0_K2_SIN_sd=sd_m_hat_mean
m_AR0_K2_SIN_CI_Inf_1=m_hat_CI[1,]
m_AR0_K2_SIN_CI_Sup_1=m_hat_CI[2,]

g11_AR0_K2_SIN_true=g11
g11_AR0_K2_SIN_esti=g11_hat_mean
g11_AR0_K2_SIN_sd=sd_g11_hat_mean
g11_AR0_K2_SIN_CI_Inf_1=g11_hat_CI[1,]
g11_AR0_K2_SIN_CI_Sup_1=g11_hat_CI[2,]

g21_AR0_K2_SIN_true=g21
g21_AR0_K2_SIN_esti=g21_hat_mean
g21_AR0_K2_SIN_sd=sd_g21_hat_mean
g21_AR0_K2_SIN_CI_Inf_1=g21_hat_CI[1,]
g21_AR0_K2_SIN_CI_Sup_1=g21_hat_CI[2,]

g12_AR0_K2_SIN_true=g12
g12_AR0_K2_SIN_esti=g12_hat_mean
g12_AR0_K2_SIN_sd=sd_g12_hat_mean
g12_AR0_K2_SIN_CI_Inf_1=g12_hat_CI[1,]
g12_AR0_K2_SIN_CI_Sup_1=g12_hat_CI[2,]

g22_AR0_K2_SIN_true=g22
g22_AR0_K2_SIN_esti=g22_hat_mean
g22_AR0_K2_SIN_sd=sd_g22_hat_mean
g22_AR0_K2_SIN_CI_Inf_1=g22_hat_CI[1,]
g22_AR0_K2_SIN_CI_Sup_1=g22_hat_CI[2,]


##############################################################################################################
##############################################################################################################
#SIMULATION 2: polynomial trend and amplitudes

w1=30*pi
w2=40*pi
lam_selec=rep(3,5)

K=2
sigma2=2

mm=0.2*x-5*x^2+5.5*x^3
g11=4*x^3-5*x^2
g12=-x+x^2+1.3*x^3
g21=-0.5-0.5*x+2.5*x^2-0.5*x^3
g22=0.5+2*x^2-3*x^3
C1=diag(cos(w1*x))
S1=diag(sin(w1*x))
C2=diag(cos(w2*x))
S2=diag(sin(w2*x))

mu=mm+C1%*%g11+S1%*%g21+C2%*%g12+S2%*%g22
plot(x,mu)

w_hat=c(w1,w2)
pordd=4
num_knots=3


set.seed(123)

coef_hat=matrix(0,ncol=SIM,nrow=((2*K+1)*(num_knots+3)))
y_hat=mu_hat=matrix(0,ncol=SIM,nrow=n)
sd_y_hat=matrix(0,ncol=SIM,nrow=n)
m_hat=sd_m_hat=matrix(0,ncol=SIM,nrow=n)
g11_hat=g21_hat=g12_hat=g22_hat=matrix(0,ncol=SIM,nrow=n)
sd_g11_hat=sd_g21_hat=sd_g12_hat=sd_g22_hat=matrix(0,ncol=SIM,nrow=n)
resid=matrix(0,ncol=SIM,nrow=n)
sigma2_hat=sigma2_hat2=rep(0,SIM)

SIM=200

for(m in 1:SIM){
  print(c(m,SIM))
  z1=rnorm(n,mean=0,sd=sqrt(sigma2))
  Y_t=mu+z1
  fit=psfit(x=x,xl=min(x)-0.0000001,xr=max(x)+0.0000001,y=Y_t,w1=w_hat,K=2,pord=pordd,
             ndx=num_knots,bdeg=3,lam=lam_selec)
  
  coef_hat[,m]=fit$theta
  y_hat[,m]=fit$f.hat
  mu_hat[,m]=fit$f.hat
  m_hat[,m]=fit$hat_m_g[,1]
  g11_hat[,m]=fit$hat_m_g[,2]
  g21_hat[,m]=fit$hat_m_g[,3]
  g12_hat[,m]=fit$hat_m_g[,4]
  g22_hat[,m]=fit$hat_m_g[,5]
  
  sigma2_hat[m]=fit$sigma2hat
  sigma2_hat2[m]=fit$sigma2hat2
  
  sd_y_hat[,m]=fit$sd_y
  sd_m_hat[,m]=fit$sd_m_g[,1]
  sd_g11_hat[,m]=fit$sd_m_g[,2]
  sd_g21_hat[,m]=fit$sd_m_g[,3]
  sd_g12_hat[,m]=fit$sd_m_g[,4]
  sd_g22_hat[,m]=fit$sd_m_g[,5]
}

#Response variable
y_hat_CI=apply(y_hat,1, quantile, probs=c(0.025,0.975))
y_hat_mean=apply(y_hat,1, mean)
sd_y_hat_mean=apply(sd_y_hat,1,mean)

mu_hat_CI=apply(mu_hat,1, quantile, probs=c(0.025,0.975))
mu_hat_mean=apply(mu_hat,1, mean)
sd_mu_hat_mean=sd_y_hat_mean

#trend and amplitudes
m_hat_CI=apply(m_hat,1, quantile, probs=c(0.025,0.975))
m_hat_mean=apply(m_hat,1, mean)
sd_m_hat_mean=apply(sd_m_hat,1,mean)
g11_hat_CI=apply(g11_hat,1, quantile, probs=c(0.025,0.975))
g11_hat_mean=apply(g11_hat,1, mean)
sd_g11_hat_mean=apply(sd_g11_hat,1, mean)
g21_hat_CI=apply(g21_hat,1, quantile, probs=c(0.025,0.975))
g21_hat_mean=apply(g21_hat,1, mean)
sd_g21_hat_mean=apply(sd_g21_hat,1, mean)
g12_hat_CI=apply(g12_hat,1, quantile, probs=c(0.025,0.975))
g12_hat_mean=apply(g12_hat,1, mean)
sd_g12_hat_mean=apply(sd_g12_hat,1, mean)
g22_hat_CI=apply(g22_hat,1, quantile, probs=c(0.025,0.975))
g22_hat_mean=apply(g22_hat,1, mean)
sd_g22_hat_mean=apply(sd_g22_hat,1, mean)

#Labels Y
Y_AR0_K2_POL=Y_t
mu_AR0_K2_POL_true=mu
mu_AR0_K2_POL_esti=mu_hat_mean
mu_AR0_K2_POL_sd=sd_mu_hat_mean
mu_AR0_K2_POL_CI_Inf_1=mu_hat_CI[1,]
mu_AR0_K2_POL_CI_Sup_1=mu_hat_CI[2,]

#Labels trend and amplitudes
m_AR0_K2_POL_true=mm
m_AR0_K2_POL_esti=m_hat_mean
m_AR0_K2_POL_sd=sd_m_hat_mean
m_AR0_K2_POL_CI_Inf_1=m_hat_CI[1,]
m_AR0_K2_POL_CI_Sup_1=m_hat_CI[2,]

g11_AR0_K2_POL_true=g11
g11_AR0_K2_POL_esti=g11_hat_mean
g11_AR0_K2_POL_sd=sd_g11_hat_mean
g11_AR0_K2_POL_CI_Inf_1=g11_hat_CI[1,]
g11_AR0_K2_POL_CI_Sup_1=g11_hat_CI[2,]

g21_AR0_K2_POL_true=g21
g21_AR0_K2_POL_esti=g21_hat_mean
g21_AR0_K2_POL_sd=sd_g21_hat_mean
g21_AR0_K2_POL_CI_Inf_1=g21_hat_CI[1,]
g21_AR0_K2_POL_CI_Sup_1=g21_hat_CI[2,]

g12_AR0_K2_POL_true=g12
g12_AR0_K2_POL_esti=g12_hat_mean
g12_AR0_K2_POL_sd=sd_g12_hat_mean
g12_AR0_K2_POL_CI_Inf_1=g12_hat_CI[1,]
g12_AR0_K2_POL_CI_Sup_1=g12_hat_CI[2,]

g22_AR0_K2_POL_true=g22
g22_AR0_K2_POL_esti=g22_hat_mean
g22_AR0_K2_POL_sd=sd_g22_hat_mean
g22_AR0_K2_POL_CI_Inf_1=g22_hat_CI[1,]
g22_AR0_K2_POL_CI_Sup_1=g22_hat_CI[2,]



pdf("SIM_K2_SIN_POLY_WN_1.pdf",width=10,height=8)

layout(matrix(c(1,1,2,2,3,3,4,4,5:12), 4, 4, byrow = TRUE),widths=c(1,1), heights=c(1,1))

par(mar = c(2,4,2,0.3))
maxx=max(mu_AR0_K2_SIN_CI_Sup_1)
minn=min(mu_AR0_K2_SIN_CI_Inf_1)
plot(x,mu_AR0_K2_SIN_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",xaxt="none",
     main='Model with sinusoidal trend and amplitudes')
lines(x,mu_AR0_K2_SIN_esti,col=1,lwd=0.5,type="l")
lines(x,mu_AR0_K2_SIN_CI_Inf_1,col=1,lty=2)
lines(x,mu_AR0_K2_SIN_CI_Sup_1,col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(mu*' and '* bar(mu)), line=2, font=2)

par(mar = c(2,4,2,0.3))
plot(x,mu_AR0_K2_POL_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",xaxt="none",yaxt="none",
     main='Model with polynomial trend and amplitudes')
lines(x,mu_AR0_K2_POL_esti,col=1,lwd=0.5,type="l")
lines(x,mu_AR0_K2_POL_CI_Inf_1,col=1,lty=2)
lines(x,mu_AR0_K2_POL_CI_Sup_1,col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(mu*' and '* bar(mu)), line=0.5, font=2)

####
p=0

maxx=2
minn=-2
par(mar = c(3,4,1,0.3))
plot(x,m_AR0_K2_SIN_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="")
lines(x,m_AR0_K2_SIN_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],m_AR0_K2_SIN_CI_Inf_1[(n*p):n*(1-p)],col=1,lwd=1,lty=2)
lines(x[(n*p):n*(1-p)],m_AR0_K2_SIN_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(m*' and '* bar(m)), line=2, font=2)

par(mar = c(3,4,1,0.3))
plot(x,m_AR0_K2_POL_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",yaxt="none")
lines(x,m_AR0_K2_POL_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],m_AR0_K2_POL_CI_Inf_1[(n*p):n*(1-p)],col=1,lwd=1,lty=2)
lines(x[(n*p):n*(1-p)],m_AR0_K2_POL_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(m*' and '* bar(m)), line=0.5, font=2)

#####

par(mar = c(3,4,1,0.3))
plot(x,g11_AR0_K2_SIN_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",xaxt="none")
lines(x,g11_AR0_K2_SIN_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],g11_AR0_K2_SIN_CI_Inf_1[(n*p):n*(1-p)],col=1,lty=2)
lines(x[(n*p):n*(1-p)],g11_AR0_K2_SIN_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[1*','*1]*' and '* bar(g)[1*','*1]), line=2, font=2)

par(mar = c(3,4,1,0.3))
plot(x,g21_AR0_K2_SIN_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",xaxt="none",yaxt="none")
lines(x,g21_AR0_K2_SIN_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],g21_AR0_K2_SIN_CI_Inf_1[(n*p):n*(1-p)],col=1,lty=2)
lines(x[(n*p):n*(1-p)],g21_AR0_K2_SIN_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[2*','*1]*' and '* bar(g)[2*','*1]), line=0.5, font=2)

#####

par(mar = c(3,4,1,0.3))
plot(x,g11_AR0_K2_POL_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",xaxt="none",yaxt="none")
lines(x,g11_AR0_K2_POL_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],g11_AR0_K2_POL_CI_Inf_1[(n*p):n*(1-p)],col=1,lty=2)
lines(x[(n*p):n*(1-p)],g11_AR0_K2_POL_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[1*','*1]*' and '* bar(g)[1*','*1]), line=0.5, font=2)

par(mar = c(3,4,1,0.3))
plot(x,g21_AR0_K2_POL_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",xaxt="none",yaxt="none")
lines(x,g21_AR0_K2_POL_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],g21_AR0_K2_POL_CI_Inf_1[(n*p):n*(1-p)],col=1,lty=2)
lines(x[(n*p):n*(1-p)],g21_AR0_K2_POL_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[2*','*1]*' and '* bar(g)[2*','*1]), line=0.5, font=2)

#####

par(mar = c(3,4,1,0.3))
plot(x,g12_AR0_K2_SIN_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="")
lines(x,g12_AR0_K2_SIN_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],g12_AR0_K2_SIN_CI_Inf_1[(n*p):n*(1-p)],col=1,lty=2)
lines(x[(n*p):n*(1-p)],g12_AR0_K2_SIN_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(g[1*','*2]*' and '* bar(g)[1*','*2]), line=2, font=2)

par(mar = c(3,4,1,0.3))
plot(x,g22_AR0_K2_SIN_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",yaxt="none")
lines(x,g22_AR0_K2_SIN_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],g22_AR0_K2_SIN_CI_Inf_1[(n*p):n*(1-p)],col=1,lty=2)
lines(x[(n*p):n*(1-p)],g22_AR0_K2_SIN_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(g[2*','*2]*' and '* bar(g)[2*','*2]), line=0.5, font=2)

####

par(mar = c(3,4,1,0.3))
plot(x,g12_AR0_K2_POL_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",yaxt="none")
lines(x,g12_AR0_K2_POL_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],g12_AR0_K2_POL_CI_Inf_1[(n*p):n*(1-p)],col=1,lty=2)
lines(x[(n*p):n*(1-p)],g12_AR0_K2_POL_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(g[1*','*2]*' and '* bar(g)[1*','*2]), line=0.5, font=2)

par(mar = c(3,4,1,0.3))
plot(x,g22_AR0_K2_POL_true,col=2,lwd=1,ylim=c(minn,maxx),type="l",xlab="",ylab="",yaxt="none")
lines(x,g22_AR0_K2_POL_esti,col=1,lwd=0.5,type="l")
lines(x[(n*p):n*(1-p)],g22_AR0_K2_POL_CI_Inf_1[(n*p):n*(1-p)],col=1,lty=2)
lines(x[(n*p):n*(1-p)],g22_AR0_K2_POL_CI_Sup_1[(n*p):n*(1-p)],col=1,lty=2)
rug(x,col='grey')
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(g[2*','*2]*' and '* bar(g)[2*','*2]), line=0.5, font=2)

dev.off()

