##########################
# AUTHOR: DARLIN SOTO    #
# DATE: OCTOBER 16, 2021 #
# EMAIL: dmsoto1@uc.cl   #
##########################

#functions

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

DFT_dar<-function(x){
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

IN_DFT_dar<-function(x){
  n=length(x)
  DFT_x=rep(0,n)
  l=1:n
  for(k in 1:n){
    DFT_x[k]=sum(x*exp(1i*k*l*2*pi/n))/N
  }
  freqs=2*pi*(1:n)/n
  result=list(DFT_x=DFT_x,freq=freqs)
  return(result)
}

###########################################################################################
##################################################################################################

#DATA
V783_sampled<-read.table("V783.txt",header=TRUE, sep = ",")
head(V783_sampled)
Y=V783_sampled[,1]
t=V783_sampled[,2]
N=length(t)

#FREQUENCIES GIVEN BY BENKO ET AL. (2014)
f0=1/0.6207001
fm=0.036058

n=15 #num harmonic 
K=1 #num hamn in the blazhko effect


############################################################################################
###########################################################################################
#FIT TO OBTAIN THE INITIAL VALUES FOR THE Levenberg-Marquardt algorithm

#DESING MATRIX
X1=matrix(0,ncol=1,nrow=N)
for(k in 1:K){
  X1=cbind(X1,cos(2*pi*k*fm*t),sin(2*pi*k*fm*t))
}
X1=X1[,-c(1)]

X2=matrix(0,ncol=1,nrow=N)
for(i in 1:n){
  X2=cbind(X2,cos(2*pi*i*f0*t),sin(2*pi*i*f0*t))
}
X2=X2[,-c(1)]

X3=matrix(0,ncol=1,nrow=N)
for(i in 1:n){
  for(j in 1:K){
    X3=cbind(X3,cos(2*pi*(i*f0+j*fm)*t),sin(2*pi*(i*f0+j*fm)*t))
  }
}
X3=X3[,-c(1)]

X4=matrix(0,ncol=1,nrow=N)
for(i in 1:n){
  for(j in 1:K){
    X4=cbind(X4,cos(2*pi*(i*f0-j*fm)*t),sin(2*pi*(i*f0-j*fm)*t))
  }
}
X4=X4[,-c(1)]

X=cbind(X1,X2,X3,X4)

fit1=lm(Y ~ X)

Y_hat=fit1$fitted.values
res=fit1$residuals

#finding parameters 
beta_hat=as.numeric(fit1$coefficients)
beta0_hat=beta_hat[1]
beta_b_hat=beta_hat[2:(2*K+1)]
beta_h_hat=beta_hat[(2*K+2):(2*K+1+2*n)]
beta_mas_hat=beta_hat[(2*K+2+2*n):(2*K+1+2*n+K*n*2)]

#calculo de phase
ind=seq(1,2*K,by=2)
phi_A=atan(beta_b_hat[ind]/beta_b_hat[-ind])
ind=seq(1,2*n,by=2)
phi_i=atan(beta_h_hat[ind]/beta_h_hat[-ind])
ind=seq(1,n*K*2,by=2)
phi_mas=atan(beta_mas_hat[ind]/beta_mas_hat[-ind])



########################################################################################################
########################################################################################################
#FIT MODEL PROPOSED BY BENKO (2018)

#INITIAL VALUES
ind=seq(1,2*K,by=2)
phi_b=atan(beta_b_hat[ind]/beta_b_hat[-ind])
ind=seq(1,2*n,by=2)
phi_i=atan(beta_h_hat[ind]/beta_h_hat[-ind])
ind=seq(1,n*K*2,by=2)
phi_mas=atan(beta_mas_hat[ind]/beta_mas_hat[-ind])

gamma_1=beta0_hat
ind=seq(1,2*K,by=2)
gamma_2=beta_b_hat[ind]/sin(phi_b)
gamma_3=phi_b
ind=seq(1,2*n,by=2)
gamma_4=beta_h_hat[ind]/sin(phi_i)
ind=seq(1,n*K*2,by=2)
gamma_5=2*beta_mas_hat[ind]/sin(phi_mas)
gamma_6=phi_mas-phi_i+pi/2
gamma_7=phi_i
gamma_8=rep(0.1,K*n)
gamma_9=rep(0.1,K*n)

#initial parameters for the Levenberg-Marquardt algorithm
gamma_int=c(gamma_1,gamma_2,gamma_3,gamma_4,gamma_5,gamma_6,gamma_7,gamma_8,gamma_9)

#Benko (2018) model
getPred <- function(gamma_int, tim){
  gamma=list()
  gamma[[1]]=gamma_int[1]
  gamma[[2]]=gamma_int[2:(1+K)]
  gamma[[3]]=gamma_int[(K+2):(2*K+1)]
  gamma[[4]]=gamma_int[(2*K+2):(2*K+1+n)]
  gamma[[5]]=gamma_int[(2*K+2+n):(2*K+1+n+K*n)]
  gamma[[6]]=gamma_int[(2*K+2+n+K*n):(2*K+1+n+2*K*n)]
  gamma[[7]]=gamma_int[(2*K+2+n+2*K*n):(2*K+1+2*n+2*K*n)]
  gamma[[8]]=gamma_int[(2*K+2+2*n+2*K*n):(2*K+1+2*n+3*K*n)]
  gamma[[9]]=gamma_int[(2*K+2+2*n+3*K*n):(2*K+1+2*n+4*K*n)]
  
  
  result=rep(0,length(tim))
  
  gamma_5_mat=matrix(gamma[[5]],ncol=K,nrow=n,byrow=TRUE)
  gamma_6_mat=matrix(gamma[[6]],ncol=K,nrow=n,byrow=TRUE)
  gamma_8_mat=matrix(gamma[[8]],ncol=K,nrow=n,byrow=TRUE)
  gamma_9_mat=matrix(gamma[[9]],ncol=K,nrow=n,byrow=TRUE)
  
  for(z in 1:length(tim)){
    term_1=gamma[[1]]
    
    term_2=sum(gamma[[2]]*sin(2*pi*(1:K)*fm*tim[z]+gamma[[3]]))
    
    term_3=0
    for(i in 1:n){
      
      aux_vect_A=matrix(0,ncol=1,nrow=K)
      for(j in 1:K){
        aux_vect_A[j]=sin(2*pi*j*fm*tim[z]+gamma_6_mat[i,j])
      }
      aux_vect_F=matrix(0,ncol=1,nrow=K)
      for(j in 1:K){
        aux_vect_F[j]=sin(2*pi*j*fm*tim[z]+gamma_9_mat[i,j])
      }
      
      aux1=gamma[[4]][i]
      aux2=gamma_5_mat[i,]%*%aux_vect_A
      aux3=gamma_8_mat[i,]%*%aux_vect_F
      aux4=sin(2*pi*i*f0*tim[z]+gamma[[7]][i]+aux3)
      term_3=term_3+(aux1+aux2)*aux4
    }
    
    result[z]=term_1+term_2+term_3 
  }
  return(result)
}

# residual function
residFun <- function(p, observed, tim) observed - getPred(p,tim)
library(minpack.lm)
#FIT
nls.out_AM_FM_2018 <- nls.lm(par=gamma_int, fn = residFun, observed = Y,tim = t)

Y_hat_AM_FM_2018=getPred(coef(nls.out_AM_FM_2018), t)
res_AM_FM_2018=residFun(coef(nls.out_AM_FM_2018), Y, t)

error_train_AM_FM_2018=mean((Y-Y_hat_AM_FM_2018)^2)

#Error
round(c(error_train_AM_FM_2018),6)
#8e-06

#time-varying trend and amplitudes
coef_hat=coef(nls.out_AM_FM_2018)
gamma_1=coef_hat[1]
gamma_2=coef_hat[2:(1+K)]
gamma_3=coef_hat[(K+2):(2*K+1)]
gamma_4=coef_hat[(2*K+2):(2*K+1+n)]
gamma_5=coef_hat[(2*K+2+n):(2*K+1+n+K*n)]
gamma_6=coef_hat[(2*K+2+n+K*n):(2*K+1+n+2*K*n)]
gamma_7=coef_hat[(2*K+2+n+2*K*n):(2*K+1+2*n+2*K*n)]
gamma_8=coef_hat[(2*K+2+2*n+2*K*n):(2*K+1+2*n+3*K*n)]
gamma_9=coef_hat[(2*K+2+2*n+3*K*n):(2*K+1+2*n+4*K*n)]

m_b=gamma_1+gamma_2*sin(2*pi*fm*t+gamma_3)
gA=gF=g1_b=g2_b=matrix(0,ncol=n,nrow=N)
for(i in 1:n){
  gA[,i]=gamma_5[i]*sin(2*pi*fm*t+gamma_6[i])
  gF[,i]=gamma_8[i]*sin(2*pi*fm*t+gamma_9[i])
  g1_b[,i]=(gamma_4[i]+gA[,i])*sin(gamma_7[i]+gF[,i])
  g2_b[,i]=(gamma_4[i]+gA[,i])*cos(gamma_7[i]+gF[,i])
}



########################################################################################################
########################################################################################################
#FIT OUR MODEL

AIC_tot=c(30,0,0.1,0.1,0,10)
K=19
K_benko=15
fe=0.6035
new_freqs=(11:14)*f0+fe
f=c((1:K_benko)*f0,new_freqs)
w_hat=2*pi*f
p=1
bdeg=3
lam_aux=c(AIC_tot[2],rep(AIC_tot[3],10),rep(AIC_tot[4],10),rep(AIC_tot[5],10),rep(AIC_tot[6],8))
knots=AIC_tot[1]


fit=psfit(x=t,xl=min(t)-0.00001,xr=max(t)+0.00001,y=Y,w1=w_hat,K=K,pord=p,ndx=knots,bdeg=bdeg,lam=lam_aux)

Y_hat_Bspline=fit$f.hat
res_Bspline=Y-Y_hat_Bspline
error_train_Bspline=mean((Y-Y_hat_Bspline)^2)

#ERROR
round(c(error_train_Bspline),6)
#1e-06

#time-varying trend and amplitudes
hat_m_g=fit$hat_m_g
sen_m_g=fit$sd_m_g


ind=seq(2,(2*K),by=2)
m_hat=hat_m_g[,1]
g1_hat=hat_m_g[,ind]
g2_hat=hat_m_g[,-c(1,ind)]  

m_sen=sen_m_g[,1]
g1_sen=sen_m_g[,ind]
g2_sen=sen_m_g[,-c(1,ind)]   


########################################################################################################
########################################################################################################
#PSD HAT OF RESIDUALS (both models)

dt_select=0.0204345
N_i=44283

#frequency
t_t0=t-min(t)
dt=dt_select
lam_l=2*pi*(1:N)/(dt*N)
lam_l=matrix(lam_l,ncol=1,nrow=length(lam_l))

#PERIODOGRAM
PER_res_AM_FM_2018=apply(lam_l,1,periodogram,x=res_AM_FM_2018,times=t_t0)
PER_res_Bspline=apply(lam_l,1,periodogram,x=res_Bspline,times=t_t0)

#SMOOTH PERIODOGRAM
s_PER_res_AM_FM_2018=ksmooth(lam_l, PER_res_AM_FM_2018, kernel ="normal", bandwidth = 7.2)$y
s_PER_res_Bspline=ksmooth(lam_l, PER_res_Bspline, kernel ="normal", bandwidth = 7.2)$y

#SPECTRAL WINDOW
W_y=apply(lam_l,1,W_fun,(t-min(t)))

#ESTIMATION PSD
aux1=DFT_dar(W_y)$DFT_x
aux2=DFT_dar(s_PER_res_AM_FM_2018)$DFT_x

aux=Re((aux2*N_i)/(aux1*2*pi))
P_res_AM_FM_2018=Re(IN_DFT_dar(aux)$DFT_x)

aux=Re((DFT_dar(s_PER_res_Bspline)$DFT_x*N)/(DFT_dar(W_y)$DFT_x*2*pi))
P_res_Bspline=Re(IN_DFT_dar(aux)$DFT_x)

ind=1:(N/2)



pdf("V783_fit_Y.pdf", width=15, height=11)

layout(matrix(1:8, 4,2, byrow = FALSE),widths=c(1,1), heights=c(1,1))

#OUT MODEL
in1=3
in2=4
in3=1
in4=1

par(mar = c(in1,in2,in3,in4))
minn_y=min(Y)-0.2
maxx_y=max(Y)+0.1
plot(t,Y,ylim=c(maxx_y,minn_y),t='l',lwd=2,col=2,ylab='',xlab='',xaxt='none')
lines(t,Y_hat_Bspline, col=1, lwd=1)
rug(t,col='grey')
title(sub='Time [BJD-2454833 d]', adj=0.5, line=0.5, font=2)
title(ylab='Brightness [mag]', line=2, font=2, cex.lab=1.3)
text(min(t)+(max(t)-min(t))/2,minn_y+0.03,cex=1.5,bquote('Fit obtained with time-varying parameters (our novel model)'))
legend("topright", legend=c('Observations','Fit'),col=c(2,"black"), lty=c(1,1), cex=1)

par(mar = c(in1,in2,in3,in4))
minn_res=min(res_AM_FM_2018,res_Bspline)
maxx_res=max(res_AM_FM_2018,res_Bspline)
plot(t,res_Bspline,ylim=c(maxx_res,minn_res),t='l',ylab='',xlab='')
rug(t,col='grey')
title(sub='Time [BJD-2454833 d]', adj=0.5, line=2, font=2)
title(ylab='Resdiduals [mag]', line=2, font=2, cex.lab=1.3)
text(min(t)+(max(t)-min(t))/2,minn_res+0.001,cex=1.3,bquote('Residual of fit obtained with time-varying parameters (our novel model)'))

maxx=max(log10(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))+0.2
minn=min(log10(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))

par(mar = c(in1,in2,in3,in4))
plot(log10(lam_l[ind]/(2*pi)),log10(P_res_Bspline[ind]),t='l',ylim=c(minn,maxx),
     xlab='',ylab='')
text(-0.15,maxx-0.2,cex=1.3,expression("log10 of estimated PSD of the residuals obtained with time-varying parameters (our novel model)"))
title(sub=bquote(paste('log10 of frequency [log10 ',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote(paste('log10 of ',hat(PSD))), line=2, font=2, cex.lab=1.3)

maxx=max(sqrt(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))+0.0008
minn=min(sqrt(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))

par(mar = c(in1,in2,in3,in4))
plot(lam_l[ind]/(2*pi),sqrt(P_res_Bspline[ind]),t='l',ylim=c(minn,maxx),xlab='',ylab='')
text(max(lam_l[ind]/(2*pi))/2,maxx-0.0003,cex=1.3,expression("Square root of estimated PSD of the residuals obtained with time-varying parameters (our novel model)"))
title(sub=bquote(paste('Frequency [',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote(sqrt(hat(PSD))), line=2, font=2, cex.lab=1.3)


#BENKO MODEL

par(mar = c(in1,in2,in3,in4))
plot(t,Y,ylim=c(maxx_y,minn_y),t='l',lwd=2,col=2,ylab='',xlab='',xaxt='none',yaxt='none')
lines(t,Y_hat_AM_FM_2018, col=1, lwd=1)
rug(t,col='grey')
title(sub='Time [BJD-2454833 d]', adj=0.5, line=0.5, font=2)
title(ylab='Brightness [mag]', line=1, font=2, cex.lab=1.3)
text(min(t)+(max(t)-min(t))/2,minn_y+0.03,cex=1.5,expression(paste('Fit obtained with time-invariant parameters (Benk',"\u00f6",' 2018)')))
#bquote('Fit obtained with Model 1'))
legend("topright", legend=c('Observations','Fit'),col=c(2,"black"), lty=c(1,1), cex=1)

par(mar = c(in1,in2,in3,in4))
plot(t,res_AM_FM_2018,ylim=c(maxx_res,minn_res),t='l',ylab='',xlab='',yaxt='none')
rug(t,col='grey')
title(sub='Time [BJD-2454833 d]', adj=0.5, line=2, font=2)
title(ylab='Resdiduals [mag]', line=1, font=2, cex.lab=1.3)
text(min(t)+(max(t)-min(t))/2,minn_res+0.001,cex=1.3,bquote('Residual of fit obtained with time-invariant parameters'))

maxx=max(log10(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))+0.2
minn=min(log10(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))

par(mar = c(in1,in2,in3,in4))
plot(log10(lam_l[ind]/(2*pi)),log10(P_res_AM_FM_2018[ind]),t='l',ylim=c(minn,maxx),
     xlab='',ylab='',yaxt='none')
text(-0.15,maxx-0.2,cex=1.3,expression("log10 of estimated PSD of the residuals obtained with time-invariant parameters"))
title(sub=bquote(paste('log10 of frequency [log10 ',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote(paste('log10 of ',hat(PSD))), line=1, font=2, cex.lab=1.3)


maxx=max(sqrt(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))+0.0008
minn=min(sqrt(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))

par(mar = c(in1,in2,in3,in4))
plot(lam_l[ind]/(2*pi),sqrt(P_res_AM_FM_2018[ind]),t='l',ylim=c(minn,maxx),xlab='',ylab='',yaxt='none')
text(max(lam_l[ind]/(2*pi))/2,maxx-0.0003,cex=1.3,expression("Square root of estimated PSD of the residuals obtained with time-invariant parameters"))
title(sub=bquote(paste('Frequency [',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote(sqrt(hat(PSD))), line=1, font=2, cex.lab=1.3)
text(x=new_freqs[1],y=0.0001,labels = bquote(paste("f"[11])))
text(x=new_freqs[1],y=0.0002,labels = bquote(paste("'")))
axis(1, at=new_freqs[1],labels = FALSE, las=1,pos=-0.0003)
text(x=new_freqs[2],y=0.0001,labels = bquote(paste("f"[12])))
text(x=new_freqs[2],y=0.0002,labels = bquote(paste("'")))
axis(1, at=new_freqs[2],labels = FALSE, las=1,pos=-0.0003)
text(x=new_freqs[3],y=0.0001,labels = bquote(paste("f"[13])))
text(x=new_freqs[3],y=0.0002,labels = bquote(paste("'")))
axis(1, at=new_freqs[3],labels = FALSE, las=1,pos=-0.0003)
text(x=new_freqs[4],y=0.0001,labels = bquote(paste("f"[14])))
text(x=new_freqs[4],y=0.0002,labels = bquote(paste("'")))
axis(1, at=new_freqs[4],labels = FALSE, las=1,pos=-0.0003)

dev.off()



pdf("V783_fit_m_g.pdf", width=14, height=10)

layout(matrix(c(rep(1,6),2:43), 8, 6, byrow = TRUE),widths=c(1,1), heights=c(1,1))
nline=0.5

index=200:1800

maxx=max(cbind(m_hat[index]+1.96*m_sen[index],m_hat[index]-1.96*m_sen[index],m_b[index]))
minn=min(cbind(m_hat[index]+1.96*m_sen[index],m_hat[index]-1.96*m_sen[index],m_b[index]))
par(mar = c(3, 3.7, 0.6, 0))
plot(t,m_hat,t='l',ylab='',xlab='',ylim=c(maxx,minn),lwd=1.5)
lines(t,m_hat+1.96*m_sen,lty=2)
lines(t,m_hat-1.96*m_sen,lty=2)
lines(t,m_b,col=2)
rug(t,col='grey')
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(m),' and ',hat(u))), line=(nline+1.5), font=2, cex.lab=1.3)


for(i in c(1,4,7,10,13)){
  par(mar = c(3, 3.7, 0.6, 0))
  maxx=max(cbind(g1_hat[index,i]+1.96*g1_sen[index,i],g1_hat[index,i]-1.96*g1_sen[index,i],g1_b[index,i]))
  minn=min(cbind(g1_hat[index,i]+1.96*g1_sen[index,i],g1_hat[index,i]-1.96*g1_sen[index,i],g1_b[index,i]))
  plot(t,g1_hat[,i],t='l',ylab='',xlab='',ylim=c(maxx,minn),xaxt='none',lwd=1.5)
  lines(t,g1_hat[,i]+1.96*g1_sen[,i],lty=2)
  lines(t,g1_hat[,i]-1.96*g1_sen[,i],lty=2)
  lines(t,g1_b[,i],col=2)
  rug(t,col='grey')
  assay <- as.character(i)
  title(sub="BJD-2450000 [d]", adj=0.5, line=0.5, font=2)
  title(ylab=bquote(paste(hat(g)[1*','*~.(assay)],' and ',hat(h)[1*','*~.(assay)])), line=(nline+1.5), font=2, cex.lab=1.3)
  
  par(mar = c(3, 3.7, 0.6, 0))
  maxx=max(cbind(g1_hat[index,i+1]+1.96*g1_sen[index,i+1],g1_hat[index,i+1]-1.96*g1_sen[index,i+1],g1_b[index,i+1]))
  minn=min(cbind(g1_hat[index,i+1]+1.96*g1_sen[index,i+1],g1_hat[index,i+1]-1.96*g1_sen[index,i+1],g1_b[index,i+1]))
  plot(t,g1_hat[,i+1],t='l',ylab='',xlab='',ylim=c(maxx,minn),xaxt='none',yaxt='none',lwd=1.5)
  lines(t,g1_hat[,i+1]+1.96*g1_sen[,i+1],lty=2)
  lines(t,g1_hat[,i+1]-1.96*g1_sen[,i+1],lty=2)
  lines(t,g1_b[,i+1],col=2)
  rug(t,col='grey')
  assay <- as.character(i+1)
  title(sub="BJD-2450000 [d]", adj=0.5, line=0.5, font=2)
  title(ylab=bquote(paste(hat(g)[1*','*~.(assay)],' and ',hat(h)[1*','*~.(assay)])), line=nline, font=2)
  
  par(mar = c(3, 3.7, 0.6, 0))
  maxx=max(cbind(g1_hat[index,i+2]+1.96*g1_sen[index,i+2],g1_hat[index,i+2]-1.96*g1_sen[index,i+2],g1_b[index,i+2]))
  minn=min(cbind(g1_hat[index,i+2]+1.96*g1_sen[index,i+2],g1_hat[index,i+2]-1.96*g1_sen[index,i+2],g1_b[index,i+2]))
  plot(t,g1_hat[,i+2],t='l',ylab='',xlab='',ylim=c(maxx,minn),xaxt='none',yaxt='none',lwd=1.5)
  lines(t,g1_hat[,i+2]+1.96*g1_sen[,i+2],lty=2)
  lines(t,g1_hat[,i+2]-1.96*g1_sen[,i+2],lty=2)
  lines(t,g1_b[,i+2],col=2)
  rug(t,col='grey')
  assay <- as.character(i+2)
  title(sub="BJD-2450000 [d]", adj=0.5, line=0.5, font=2)
  title(ylab=bquote(paste(hat(g)[1*','*~.(assay)],' and ',hat(h)[1*','*~.(assay)])), line=nline, font=2, cex.lab=1.3)
  

  par(mar = c(3, 3.7, 0.6, 0))
  maxx=max(cbind(g2_hat[index,i]+1.96*g2_sen[index,i],g2_hat[index,i]-1.96*g2_sen[index,i],g2_b[index,i]))
  minn=min(cbind(g2_hat[index,i]+1.96*g2_sen[index,i],g2_hat[index,i]-1.96*g2_sen[index,i],g2_b[index,i]))
  plot(t,g2_hat[,i],t='l',ylab='',xlab='',ylim=c(maxx,minn),xaxt='none',yaxt='none',lwd=1.5)
  lines(t,g2_hat[,i]+1.96*g2_sen[,i],lty=2)
  lines(t,g2_hat[,i]-1.96*g2_sen[,i],lty=2)
  lines(t,g2_b[,i],col=2)
  rug(t,col='grey')
  assay <- as.character(i)
  title(sub="BJD-2450000 [d]", adj=0.5, line=0.5, font=2)
  title(ylab=bquote(paste(hat(g)[2*','*~.(assay)],' and ',hat(h)[2*','*~.(assay)])), line=nline, font=2, cex.lab=1.3)
  
  par(mar = c(3, 3.7, 0.6, 0))
  maxx=max(cbind(g2_hat[index,i+1]+1.96*g2_sen[index,i+1],g2_hat[index,i+1]-1.96*g2_sen[index,i+1],g2_b[index,i+1]))
  minn=min(cbind(g2_hat[index,i+1]+1.96*g2_sen[index,i+1],g2_hat[index,i+1]-1.96*g2_sen[index,i+1],g2_b[index,i+1]))
  plot(t,g2_hat[,i+1],t='l',ylab='',xlab='',ylim=c(maxx,minn),xaxt='none',yaxt='none',lwd=1.5)
  lines(t,g2_hat[,i+1]+1.96*g2_sen[,i+1],lty=2)
  lines(t,g2_hat[,i+1]-1.96*g2_sen[,i+1],lty=2)
  lines(t,g2_b[,i+1],col=2)
  rug(t,col='grey')
  assay <- as.character(i+1)
  title(sub="BJD-2450000 [d]", adj=0.5, line=0.5, font=2)
  title(ylab=bquote(paste(hat(g)[2*','*~.(assay)],' and ',hat(h)[2*','*~.(assay)])), line=nline, font=2, cex.lab=1.3)
  
  par(mar = c(3, 3.7, 0.6, 0))
  maxx=max(cbind(g2_hat[index,i+2]+1.96*g2_sen[index,i+2],g2_hat[index,i+2]-1.96*g2_sen[index,i+2],g2_b[index,i+2]))
  minn=min(cbind(g2_hat[index,i+2]+1.96*g2_sen[index,i+2],g2_hat[index,i+2]-1.96*g2_sen[index,i+2],g2_b[index,i+2]))
  plot(t,g2_hat[,i+2],t='l',ylab='',xlab='',ylim=c(maxx,minn),xaxt='none',yaxt='none',lwd=1.5)
  lines(t,g2_hat[,i+2]+1.96*g2_sen[,i+2],lty=2)
  lines(t,g2_hat[,i+2]-1.96*g2_sen[,i+2],lty=2)
  lines(t,g2_b[,i+2],col=2)
  rug(t,col='grey')
  assay <- as.character(i+2)
  title(sub="BJD-2450000 [d]", adj=0.5, line=0.5, font=2)
  title(ylab=bquote(paste(hat(g)[2*','*~.(assay)],' and ',hat(h)[2*','*~.(assay)])), line=nline, font=2, cex.lab=1.3)

}

i=11

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g1_hat[index,i]+1.96*g1_sen[index,i],g1_hat[index,i]-1.96*g1_sen[index,i]))
minn=min(cbind(g1_hat[index,i]+1.96*g1_sen[index,i],g1_hat[index,i]-1.96*g1_sen[index,i]))
plot(t,g1_hat[,i],t='l',ylab='',xlab='',ylim=c(maxx,minn),lwd=1.5)
lines(t,g1_hat[,i]+1.96*g1_sen[,i],lty=2)
lines(t,g1_hat[,i]-1.96*g1_sen[,i],lty=2)
rug(t,col='grey')
assay <- as.character(i)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[1*','*~.(assay)])), line=(nline+1.5), font=2, cex.lab=1.3)


par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g1_hat[index,i+1]+1.96*g1_sen[index,i+1],g1_hat[index,i+1]-1.96*g1_sen[index,i+1]))
minn=min(cbind(g1_hat[index,i+1]+1.96*g1_sen[index,i+1],g1_hat[index,i+1]-1.96*g1_sen[index,i+1]))
plot(t,g1_hat[,i+1],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g1_hat[,i+1]+1.96*g1_sen[,i+1],lty=2)
lines(t,g1_hat[,i+1]-1.96*g1_sen[,i+1],lty=2)
rug(t,col='grey')
assay <- as.character(i+1)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[1*','*~.(assay)])), line=(nline), font=2, cex.lab=1.3)

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g1_hat[index,i+2]+1.96*g1_sen[index,i+2],g1_hat[index,i+2]-1.96*g1_sen[index,i+2]))
minn=min(cbind(g1_hat[index,i+2]+1.96*g1_sen[index,i+2],g1_hat[index,i+2]-1.96*g1_sen[index,i+2]))
plot(t,g1_hat[,i+2],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g1_hat[,i+2]+1.96*g1_sen[,i+2],lty=2)
lines(t,g1_hat[,i+2]-1.96*g1_sen[,i+2],lty=2)
rug(t,col='grey')
assay <- as.character(i+2)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[1*','*~.(assay)])), line=(nline), font=2, cex.lab=1.3)

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g2_hat[index,i]+1.96*g2_sen[index,i],g2_hat[index,i]-1.96*g2_sen[index,i]))
minn=min(cbind(g2_hat[index,i]+1.96*g2_sen[index,i],g2_hat[index,i]-1.96*g2_sen[index,i]))
plot(t,g2_hat[,i],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g2_hat[,i]+1.96*g2_sen[,i],lty=2)
lines(t,g2_hat[,i]-1.96*g2_sen[,i],lty=2)
rug(t,col='grey')
assay <- as.character(i)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[2*','*~.(assay)])), line=(nline), font=2, cex.lab=1.3)

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g2_hat[index,i+1]+1.96*g2_sen[index,i+1],g2_hat[index,i+1]-1.96*g2_sen[index,i+1]))
minn=min(cbind(g2_hat[index,i+1]+1.96*g2_sen[index,i+1],g2_hat[index,i+1]-1.96*g2_sen[index,i+1]))
plot(t,g2_hat[,i+1],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g2_hat[,i+1]+1.96*g2_sen[,i+1],lty=2)
lines(t,g2_hat[,i+1]-1.96*g2_sen[,i+1],lty=2)
rug(t,col='grey')
assay <- as.character(i+1)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[2*','*~.(assay)])), line=(nline), font=2, cex.lab=1.3)

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g2_hat[index,i+2]+1.96*g2_sen[index,i+2],g2_hat[index,i+2]-1.96*g2_sen[index,i+2]))
minn=min(cbind(g2_hat[index,i+2]+1.96*g2_sen[index,i+2],g2_hat[index,i+2]-1.96*g2_sen[index,i+2]))
plot(t,g2_hat[,i+2],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g2_hat[,i+2]+1.96*g2_sen[,i+2],lty=2)
lines(t,g2_hat[,i+2]-1.96*g2_sen[,i+2],lty=2)
rug(t,col='grey')
assay <- as.character(i+2)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[2*','*~.(assay)])), line=(nline), font=2, cex.lab=1.3)

i=14

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g1_hat[index,i]+1.96*g1_sen[index,i],g1_hat[index,i]-1.96*g1_sen[index,i]))
minn=min(cbind(g1_hat[index,i]+1.96*g1_sen[index,i],g1_hat[index,i]-1.96*g1_sen[index,i]))
plot(t,g1_hat[,i],t='l',ylab='',xlab='',ylim=c(maxx,minn),lwd=1.5)
lines(t,g1_hat[,i]+1.96*g1_sen[,i],lty=2)
lines(t,g1_hat[,i]-1.96*g1_sen[,i],lty=2)
rug(t,col='grey')
assay <- as.character(i)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[1*','*~.(assay)])), line=(nline+1.5), font=2, cex.lab=1.3)

par(mar = c(3, 3.7, 0.6, 0))
plot.new()
par(mar = c(3, 3.7, 0.6, 0))
plot.new()

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g2_hat[index,i]+1.96*g2_sen[index,i],g2_hat[index,i]-1.96*g2_sen[index,i]))
minn=min(cbind(g2_hat[index,i]+1.96*g2_sen[index,i],g2_hat[index,i]-1.96*g2_sen[index,i]))
plot(t,g2_hat[,i],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g2_hat[,i]+1.96*g2_sen[,i],lty=2)
lines(t,g2_hat[,i]-1.96*g2_sen[,i],lty=2)
rug(t,col='grey')
assay <- as.character(i)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[2*','*~.(assay)])), line=(nline), font=2, cex.lab=1.3)

par(mar = c(3, 3.7, 0.6, 0))
plot.new()
par(mar = c(3, 3.7, 0.6, 0))
plot.new()

dev.off()

