##########################
# AUTHOR: DARLIN SOTO    #
# DATE: OCTOBER 16, 2021 #
# EMAIL: dmsoto1@uc.cl   #
##########################

library(pracma)

#FUNCTIONS
spectral_density <- function(ar = numeric(), ma = numeric(), sd = 1, lambda = NULL) {
  p <- length(ar)
  q <- length(ma)
  phi <- c(1, -ar)
  theta <- c(1, +ma)
  if (is.null(lambda)) {
    lambda <- seq(0, pi, 0.01)
  }
  n <- length(lambda)
  aux <- c()
  for (k in 1:n) {
    aux[k] <- (sum((theta * exp(-1i * lambda[k] * c(0:q)))) * sum((theta * exp(+1i * lambda[k] * c(0:q))))) / (sum((phi * exp(-1i * lambda[k] * c(0:p)))) * sum((phi * exp(+1i * lambda[k] * c(0:p)))))
  }
  sigma <- sd
  aux <- sigma^2 * Re(aux) / (2 * pi)
  return(list(PSD = aux, lambda = lambda))
}

W_fun<-function(freq,time){
  aux=sum(exp(1i*freq*time))
  aux_conv=sum(exp(-1i*freq*time))
  Re(aux*aux_conv)
}

periodogram<-function(lam,x,times){
  A = sum(x*cos(times*lam))
  B = sum(x*sin(times*lam))
  A^2+B^2
}

DFT<-function(x){
  n=length(x)
  DFT_x=rep(0,n)
  l=1:n
  for(k in 1:n){
    DFT_x[k]=sum(x*exp(-1i*k*l*2*pi/n))
  }
  freqs=2*pi*(1:n)/n
  result=list(DFT_x=DFT_x,freq=freqs)
  return(result)
}

IN_DFT<-function(x){
  n=length(x)
  DFT_x=rep(0,n)
  l=1:n
  for(k in 1:n){
    DFT_x[k]=sum(x*exp(1i*k*l*2*pi/n))/n
  }
  freqs=2*pi*(1:n)/n
  result=list(DFT_x=DFT_x,freq=freqs)
  return(result)
}


#########################################################################################
#########################################################################################

#parameters AR(2)
N=500
dt=0.33
times_x=seq(1,by=dt,length=N)
phi1=1.318
phi2=-0.634
sigma2=289.2

#sample
Num_block=50
Num_obs_block=N/Num_block
ind_block=sort(rep(1:Num_block,Num_obs_block))
sort(sample(1:Num_block,30,replace = FALSE))
#2  4  5  6  7  8 10 17 19 20 21 22 24 25 26 28 29 30 31 32 34 35 37 39 40 41 42 44 47 49

index_y=which(ind_block==2 | ind_block==4 | ind_block==5 | ind_block==6 | ind_block==7 | ind_block==8 
             | ind_block==10 | ind_block==17 | ind_block==19 | ind_block==20 | ind_block==21 
             | ind_block==22 | ind_block==24 | ind_block==25 | ind_block==26 | ind_block==28
             | ind_block==29 | ind_block==30 | ind_block==31 | ind_block==32 | ind_block==34
             | ind_block==35 | ind_block==37 | ind_block==39 | ind_block==40 | ind_block==41
             | ind_block==42 | ind_block==44 | ind_block==47 | ind_block==50)

N_y=max(index_y)
times_y=times_x[index_y]

#frequency
lam_l_x=2*pi*(1:N)/(dt*N)
N_l_x=length(lam_l_x)
lam_l_x_vect=matrix(lam_l_x,ncol=1,nrow=N_l_x)
lam_l_y=2*pi*(1:N_y)/(dt*N_y)
N_l_y=length(lam_l_y)
lam_l_y_vect=matrix(lam_l_y,ncol=1,nrow=N_l_y)


SIM=200
PER_y=PER_x=matrix(0,ncol=SIM,nrow=N_l_y)

set.seed(1)

for(m in 1:SIM){
  z=rnorm(N,mean=0,sd=sqrt(sigma2))
  x=z

  for(i in 3:N){
    x[i]=phi1*x[i-1]+phi2*x[i-2]+z[i]
  }  
  
  PER_x[,m]=apply(lam_l_y_vect,1,periodogram,x=x,times=times_x)
  
  y=x[index_y]
  PER_y[,m]=apply(lam_l_y_vect,1,periodogram,x=y,times=times_y)
}

PER_mean_y=apply(PER_y,1,mean)
PER_mean_x=apply(PER_x,1,mean)

#TRUE SPECTRAL DENSITY
aux=spectral_density(ar = c(phi1,phi2), ma = numeric(), sd = sqrt(sigma2), lambda = lam_l_x_vect*dt)
PSD=aux$PSD

plot(lam_l_x_vect[1:250],PER_mean_x[1:250]/500,t='l',col='gray')
lines(lam_l_y_vect[1:250],PER_mean_y[1:250]/300,t='l',col='azure4')
lines(lam_l_x_vect[1:250],2*pi*PSD[1:250],col=1,lty=2)

#ESTIMATED SPECTRAL DENSITY
W_y=apply(lam_l_y_vect,1,W_fun,times_y)
aux=Re((DFT(PER_mean_y)$DFT_x*N_y)/(DFT(W_y)$DFT_x*2*pi))
P_y_hat=Re(IN_DFT(aux)$DFT_x)

s_P_y_hat=ksmooth(lam_l_y_vect, P_y_hat, kernel ="normal", bandwidth = 0.3)$y

ind=1:(N_y/2)
f=(1:N)/(dt*N)


pdf("SIM_PSD_AR2.pdf",width=12,height=4)

layout(matrix(c(1,2,3), 1, 3, byrow = TRUE),widths=c(1,1), heights=c(1,1))

maxx=max(P_y_hat[ind])
minn=min(P_y_hat[ind])

par(mar=c(3.3, 3.3, 1, 1))
plot(f[ind],PER_mean_x[ind]/(2*pi*500),t='l',col='darkgreen',ylim=c(minn,maxx),ylab='',xlab='')
lines(f[ind],PER_mean_y[ind]/(2*pi*300),t='l',col='darkblue')
lines(f[ind],PSD[ind],col=2,lwd=1)
abline(h=0)
title(sub="Frequency", adj=0.5, line=2, font=2)
title(ylab='Spectral density', line=2, font=2)
legend("topright", legend=c(as.expression(bquote(bar(I)[epsilon]/(2*pi*N[y]))),
                            as.expression(bquote(bar(I)[e]/(2*pi*N[x]))),
                            as.expression(bquote(P[epsilon]))),col=c('darkgreen','darkblue',2), lwd=2,lty=c(1,1), cex=1)
text(max(f[ind])/2,maxx,cex=1.3,bquote('Periodograms'))

par(mar=c(3.3, 3.3, 1, 1))
plot(f[ind],P_y_hat[ind],col='blue',lwd=1,t='l',ylim=c(minn,maxx),ylab='',xlab='',yaxt="none")
lines(f[ind],PSD[ind],col=2,lwd=1)
abline(h=0)
title(sub="Frequency", adj=0.5, line=2, font=2)
title(ylab='Spectral density', line=0.5, font=2)
legend("topright", legend=c(as.expression(bquote(hat(P)[e])),
                            as.expression(bquote(P[epsilon]))),col=c('blue', 2),lwd=2, lty=c(1,1), cex=1)
text(max(f[ind])/2,maxx,cex=1.3,bquote('Fit obtained without smoothing'))

par(mar=c(3.3, 3.3, 1, 1))
plot(f[ind],s_P_y_hat[ind],col=2,lwd=1,t='l',ylim=c(minn,maxx),ylab='',xlab='',yaxt="none")
lines(f[ind],PSD[ind],col='blue',lwd=1)
abline(h=0)
title(sub="Frequency", adj=0.5, line=2, font=2)
title(ylab='Spectral density', line=0.5, font=2)
legend("topright", legend=c(as.expression(bquote(tilde(P)[e])),
                            as.expression(bquote(P[epsilon]))),col=c('blue', 2), lwd=2,lty=c(1,1), cex=1)
text(max(f[ind])/2,maxx,cex=1.3,bquote('Fit obtained with bandwidth=0.3'))

dev.off()
