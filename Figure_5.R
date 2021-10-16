rm(list=ls())

setwd("~/Doctorado/Tesis/Modulation models/Paper ApJ 2021/R code Apj 2021 WN errors")

library(lomb)
library(splines)
source("~/Doctorado/Tesis/Modulation models/Paper ApJ 2021/R code Apj 2021 WN errors/B_spline_functions/psfit3.R")
source("~/Doctorado/Tesis/Modulation models/Paper ApJ 2021/R code Apj 2021 WN errors/B_spline_functions/psfit3_AIC.R")
source("~/Doctorado/Tesis/Modulation models/Paper ApJ 2021/R code Apj 2021 WN errors/B_spline_functions/bspline.R")
source("~/Doctorado/Tesis/Modulation models/Paper ApJ 2021/R code Apj 2021 WN errors/B_spline_functions/bspredict.R")

memory.limit(size=30000)
options(digits=11)

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



star<-read.delim('V783_edit.txt', header =FALSE, sep = ";")

tt=star[,3]
yy=star[,8]

F_n=24.46842


#per=lsp(yy,times=tt,from=0,to=2*F_n,type='frequency',ofac=1)
#per$peak.at
#f0[1]=per$peak.at[1]

#pdf("PER_WHOLE_V783.pdf", width=10, height=7)
#plot(per$scanned,per$power,col=1,t='l',xlab='Frequency',ylab='Periodogram')
#abline(v=F_n,col=2)
#dev.off()



#N=3000
#f_min=0.01
#f_max=200
#v=seq(f_min,f_max,length=N)

#lambda=2*pi*v
#lambda=matrix(lambda,ncol=1,nrow=N)

#SW=apply(lambda,1,W_fun,time=tt)

#pdf("SW_WHOLE_V783.pdf", width=10, height=7)
#plot(v,SW/length(tt),t='l',ylab='Spectral window',xlab='Frequency')
#dev.off()




head(star)
dim(star)
summary(star)
nrow(star)

#SELECCIONAMOS DATOS CON ESTRUCTURA t_i=i*delta
aux_round=7
data_t=round(star[,3],aux_round)
data_t
sum(table(data_t)==1)
N_T=length(star[,3])

dt=diff(data_t)
data_t=cbind(data_t,c(1,dt))
data_t
data_t=round(data_t,aux_round)

#select dt
dt_round=round(dt,aux_round)
dt_data <- data.frame(dt = dt_round)
new_df<-aggregate(dt_data$dt, dt_data, length)
new_df
ind=which.max(new_df[,2])
ind
new_df[ind,]#0.0204345
dt_select=new_df[ind,1]
data_t=cbind(data_t,rep(dt_select,nrow(data_t)))
data_t=round(data_t,aux_round)

#dt_select
data_t=cbind(data_t,rep(0,nrow(data_t)))
ind_t=which(data_t[,2]==dt_select)
ind_t
data_t[ind_t,4]=1
data_t[ind_t-1,4]=1
data_t
i_index=data_t[,1]/data_t[,3]
data_t=cbind(data_t,i_index)
data_t=round(data_t,aux_round)
data_t
t_new1=data_t[which(data_t[,4]==1),1]
i_ind1=data_t[which(data_t[,4]==1),5]
i_ind1
plot(i_ind1,t='l')

dec <- data_t[,5] - floor(data_t[,5])
dec
plot(dec,t='l')
plot(dec[which(data_t[,4]==1)],t='l',col=2)
data_t=cbind(data_t,dec)
data_t=round(data_t,aux_round)
head(data_t)

aux=dec[which(data_t[,4]==1)]
plot(aux,t='l')
head(data_t[which(data_t[,4]==1),])
plot(aux[3500:5600],t='l')
min(aux[3500:5600])
max(aux[3500:5600])
length(aux[3500:5600])
index_fil=3500:5600
#2101
#filtro despues de 4

#creamos nuevo tiempo desplazado para ener enteros
aux=min(aux[index_fil])
t0=aux*dt_select
t0#0.0035986993605
t_new2=data_t[,1]-t0
data_t=cbind(data_t,t_new2)
data_t=round(data_t,aux_round)
head(data_t)
i_index2=data_t[,7]/data_t[,3]
data_t=cbind(data_t,i_index2)
data_t=round(data_t,aux_round)
head(data_t)


#Seleccionamos observaciones con i_index2==entero
data_t=cbind(data_t,star[,8])
data_t=round(data_t,aux_round)
head(data_t)
aux=which(data_t[,4]==1)
length(aux)

data1=data_t[aux,]
nrow(data1)
data2=data1[index_fil,]
nrow(data2)
#2101 obs

star_new=star[aux,]
star_new=star_new[index_fil,]
dim(star_new)

star_new=star_new[,-1]
write.table(star_new, file = "V783_sample2.txt", append = FALSE, quote = TRUE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)



t<-data2[,1]
Y<-data2[,9]
N=length(t)
N
N_i=round(max(data2[,8]),1)
N_i
plot(t,Y,t='l')

round(max(t)-min(t),0)#77
round(max(t),2)#904.9
round(min(t),1)#827.4

i_ind=round(data2[,8],1)
dt=dt_select
t_new=t0+i_ind*dt
plot(t[1:10],t='l',lwd=2)
lines(t_new[1:10],col=2)

t0_new2=t0+40491*dt
t_new_2=t0_new2+(i_ind-40491)*dt

t_new[1:10]
t_new_2[1:10]

write.table(i_ind, file = "i_ind_sample2.txt", append = FALSE, quote = TRUE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(t_new, file = "t_sample2.txt", append = FALSE, quote = TRUE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)


del=3
lam_1=2*pi*(1:N_i)/(dt*N_i)
lam_1=matrix(lam_1,ncol=1,nrow=length(lam_1))
lam_2=2*pi*((1-del):(N_i-del))/(dt*N_i)
lam_2=matrix(lam_2,ncol=1,nrow=length(lam_2))
W_1=apply(lam_1,1,W_fun,t_new)
W_2=apply(lam_2,1,W_fun,t_new)
c(sum(W_1),sum(W_2))

plot(t,Y,t='l',lwd=2)
lines(t_new,Y,col=2)

t=t_new
V783_sample2_R=cbind(t,Y)
dim(V783_sample2_R)

write.table(V783_sample2_R, file = "V783_sample2_R.txt", append = FALSE, quote = TRUE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)


per=lsp(Y,from = 0, to = 14,type='frequency',ofac=2)
per$peak.at
PERIOD=per$power
freq=per$scanned
plot(freq,PERIOD,t='l',ylim=c(0,10),col=2,lty=2)



f0=1/0.6207001
fm=0.036058

n=15 #num harmonic 
K=1 #num hamn in the blazhko effect

round(f0,6)#1.611084
round(fm,6)#0.036058

################ AJUSTE POR OLS PARA OBTENER PARAMETROS INICIALES ##########################

#DESING MATRIX
X1=matrix(0,ncol=1,nrow=N)
for(k in 1:K){
  X1=cbind(X1,cos(2*pi*k*fm*t),sin(2*pi*k*fm*t))
}
X1=X1[,-c(1)]
head(X1)
dim(X1)

X2=matrix(0,ncol=1,nrow=N)
for(i in 1:n){
  X2=cbind(X2,cos(2*pi*i*f0*t),sin(2*pi*i*f0*t))
}
X2=X2[,-c(1)]
head(X2)
dim(X2)

X3=matrix(0,ncol=1,nrow=N)
for(i in 1:n){
  for(j in 1:K){
    X3=cbind(X3,cos(2*pi*(i*f0+j*fm)*t),sin(2*pi*(i*f0+j*fm)*t))
  }
}
X3=X3[,-c(1)]
head(X3)
dim(X3)

X4=matrix(0,ncol=1,nrow=N)
for(i in 1:n){
  for(j in 1:K){
    X4=cbind(X4,cos(2*pi*(i*f0-j*fm)*t),sin(2*pi*(i*f0-j*fm)*t))
  }
}
X4=X4[,-c(1)]
head(X4)
dim(X4)

X=cbind(X1,X2,X3,X4)
dim(X)

fit1=lm(Y ~ X)
summary(fit1)

Y_hat=fit1$fitted.values
res=fit1$residuals

plot(t,Y,ylim=c(max(Y),min(Y)),t='l')
lines(t,Y_hat,col=2)

plot(t,res,t='l')


#finding parameters 

beta_hat=as.numeric(fit1$coefficients)
length(beta_hat)

beta0_hat=beta_hat[1]

beta_b_hat=beta_hat[2:(2*K+1)]
length(beta_b_hat)
dim(X1)

beta_h_hat=beta_hat[(2*K+2):(2*K+1+2*n)]
length(beta_h_hat)
dim(X2)

beta_mas_hat=beta_hat[(2*K+2+2*n):(2*K+1+2*n+K*n*2)]
length(beta_mas_hat)
dim(X3)

#calculo de phase
ind=seq(1,2*K,by=2)
phi_A=atan(beta_b_hat[ind]/beta_b_hat[-ind])

ind=seq(1,2*n,by=2)
phi_i=atan(beta_h_hat[ind]/beta_h_hat[-ind])

ind=seq(1,n*K*2,by=2)
phi_mas=atan(beta_mas_hat[ind]/beta_mas_hat[-ind])



########################################################################################################
########################################################################################################
########################################################################################################
############################################   BENKO 2018   ############################################
########################################################################################################
########################################################################################################
########################################################################################################


#total of initial values for BENKO 2011

alpha_1=beta0_hat
ind=seq(1,2*K,by=2)
alpha_2=beta_b_hat[ind]/sin(phi_A)
alpha_3=phi_A
ind=seq(1,2*n,by=2)
alpha_4=beta_h_hat[ind]/sin(phi_i)
ind=seq(1,n*K*2,by=2)
alpha_5=2*beta_mas_hat[ind]/sin(phi_mas)
alpha_6=phi_i
alpha_7=rep(0,K)
alpha_8=rep(0,K)


#calculo de phase
ind=seq(1,2*K,by=2)
phi_b=atan(beta_b_hat[ind]/beta_b_hat[-ind])

ind=seq(1,2*n,by=2)
phi_i=atan(beta_h_hat[ind]/beta_h_hat[-ind])

ind=seq(1,n*K*2,by=2)
phi_mas=atan(beta_mas_hat[ind]/beta_mas_hat[-ind])

#total of initial values for BENKO 2018

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


library(minpack.lm)


# MODELO 1: BENKO 2018 AM AND FM

gamma_int=c(gamma_1,gamma_2,gamma_3,gamma_4,gamma_5,gamma_6,gamma_7,gamma_8,gamma_9)

## model based on a list of parameters
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

check <- getPred(gamma_int,t) 
plot(t,Y,ylim=c(max(Y),min(Y)),t='l')
lines(t,check,col=2)

plot(t,check,t='l')

## residual function
residFun <- function(p, observed, tim) observed - getPred(p,tim)
residFun(gamma_int, Y, t)

#start_time <- Sys.time()
#nls.out_AM_FM_2018 <- nls.lm(par=gamma_int, fn = residFun, observed = Y,
#                             tim = t)
#end_time <- Sys.time()
#end_time - start_time


#gamma_int=rep(0.1,93)
## perform fit
nls.out_AM_FM_2018 <- nls.lm(par=gamma_int, fn = residFun, observed = Y,
                             tim = t)
nls.out_AM_FM_2018$message
Y_hat_AM_FM_2018=getPred(coef(nls.out_AM_FM_2018), t)
write.table(Y_hat_AM_FM_2018, file = "V783_Y_hat_AM_FM_2018_SAMPLE2.txt", append = FALSE, quote = TRUE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
Y_hat_AM_FM_2018<-as.matrix(read.delim('V783_Y_hat_AM_FM_2018_SAMPLE2.txt', header =FALSE, sep = " "))

plot(t,Y,ylim=c(max(Y),min(Y)),t='l',lwd=2)
lines(t,Y_hat_AM_FM_2018, col=2, lwd=1)

res_AM_FM_2018=residFun(coef(nls.out_AM_FM_2018), Y, t)
write.table(res_AM_FM_2018, file = "V783_res_AM_FM_2018_SAMPLE2.txt", append = FALSE, quote = TRUE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
res_AM_FM_2018<-as.matrix(read.delim('V783_res_AM_FM_2018_SAMPLE2.txt', header =FALSE, sep = " "))



F_n=(1/(2*dt))
per=lsp(res_AM_FM_2018,times=t,from=0,to=F_n,type='frequency',ofac=1)
period=per$power
freq=per$scanned

ind=which.max(period)
peak1=c(freq[ind],period[ind])
peak1
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak2=c(freq[ind],period[ind])
peak2
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak3=c(freq[ind],period[ind])
peak3
period[1:(ind+10)]=0
ind=which.max(period)
peak4=c(freq[ind],period[ind])
peak4
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak5=c(freq[ind],period[ind])
peak5
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak6=c(freq[ind],period[ind])
peak6
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak7=c(freq[ind],period[ind])
peak7
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak8=c(freq[ind],period[ind])
peak8


per=lsp(res_AM_FM_2018,times=t,from=25.7,to=2*F_n,type='frequency',ofac=1)
period=per$power
freq=per$scanned

ind=which.max(period)
peak1_n=c(freq[ind],period[ind])
peak1_n
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak2_n=c(freq[ind],period[ind])
peak2_n
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak3_n=c(freq[ind],period[ind])
peak3_n
period[1:(ind+10)]=0
ind=which.max(period)
peak4_n=c(freq[ind],period[ind])
peak4_n
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak5_n=c(freq[ind],period[ind])
peak5_n
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak6_n=c(freq[ind],period[ind])
peak6_n
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak7_n=c(freq[ind],period[ind])
peak7_n
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak8_n=c(freq[ind],period[ind])
peak8_n
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak9_n=c(freq[ind],period[ind])
peak9_n
period[(ind-10):(ind+10)]=0
ind=which.max(period)
peak10_n=c(freq[ind],period[ind])
peak10_n

n=1:14
F_n=24.46842
F_NEW=2*F_n-(30-n)*f0


per=lsp(res_AM_FM_2018,times=t,from=0,to=2*F_n,type='frequency',ofac=1)

siz=0.9

pdf("PER1_RESIDUAL_BENKO.pdf", width=15, height=9)

par(mar = c(4.3, 4, 0.5, 0.5))
plot(per$scanned,per$power,col=1,t='l',xlab='Frequency',ylab='Periodogram')
abline(v=F_n,col=2)

#freqs left
text(x=peak1[1],y=0.004+peak1[2],labels = bquote(paste(f[14],'=2',f[N],'-16',f[0])),
     cex=siz)
text(x=peak1[1]-1.3,y=0.004+peak1[2]+0.0008,labels = bquote(paste("'")),
     cex=siz)
text(x=peak2[1],y=0.004+peak2[2],labels = bquote(paste(f[13],'=2',f[N],'-17',f[0])),
     cex=siz)
text(x=peak2[1]-1.3,y=0.004+peak2[2]+0.0008,labels = bquote(paste("'")),
     cex=siz)
text(x=peak4[1],y=0.004+peak4[2],labels = bquote(paste(f[12],'=2',f[N],'-18',f[0])),
     cex=siz)
text(x=peak4[1]-1.3,y=0.004+peak4[2]+0.0008,labels = bquote(paste("'")),
     cex=siz)
text(x=peak5[1],y=0.004+peak5[2],labels = bquote(paste(f[10],'=2',f[N],'-20',f[0])),
     cex=siz)
text(x=peak5[1]-1.3,y=0.004+peak5[2]+0.0008,labels = bquote(paste("'")),
     cex=siz)
text(x=peak6[1],y=0.001+peak6[2],labels = bquote(paste(f[11],'=2',f[N],'-19',f[0])),
     cex=siz)
text(x=peak6[1]-1.3,y=0.0009+peak6[2]+0.0008,labels = bquote(paste("'")),
     cex=siz)
text(x=peak7[1],y=0.004+peak7[2],labels = bquote(paste(f[9],'=2',f[N],'-21',f[0])),
     cex=siz)
text(x=peak7[1]-1.2,y=0.004+peak7[2]+0.0008,labels = bquote(paste("'")),
     cex=siz)
text(x=peak8[1],y=0.004+peak8[2],labels = bquote(paste(f[8],'=2',f[N],'-22',f[0])),
     cex=siz)
text(x=peak8[1]-1.2,y=0.004+peak8[2]+0.0008,labels = bquote(paste("'")),
     cex=siz)

#freqs right
text(x=peak1_n[1],y=0.004+peak1_n[2],labels = bquote(paste('16',f[0])),
     cex=siz)
text(x=peak4_n[1],y=0.004+peak4_n[2],labels = bquote(paste('17',f[0])),
     cex=siz)
text(x=peak6_n[1],y=0.004+peak6_n[2],labels = bquote(paste('18',f[0])),
     cex=siz)
text(x=peak7_n[1],y=0.004+peak7_n[2],labels = bquote(paste('20',f[0])),
     cex=siz)
text(x=peak8_n[1],y=0.004+peak8_n[2],labels = bquote(paste('19',f[0])),
     cex=siz)
text(x=peak9_n[1],y=0.004+peak9_n[2],labels = bquote(paste('21',f[0])),
     cex=siz)
text(x=peak10_n[1],y=0.004+peak10_n[2],labels = bquote(paste('22',f[0])),
     cex=siz)

dev.off()









error_train_AM_FM_2018=mean((Y-Y_hat_AM_FM_2018)^2)

plot(t,res_AM_FM_2018,t='l')

round(c(error_train_AM_FM_2018),6)
#8e-06

#g^A, g^F
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

#frequency
t_t0=t-min(t)
t_int=t_t0/dt

dt=dt_select
lam_l=2*pi*(1:N)/(N)
lam_l=matrix(lam_l,ncol=1,nrow=length(lam_l))
ind=1:(N/2)


PER_Y=apply(lam_l,1,periodogram,x=res_AM_FM_2018,times=t_int)
plot(lam_l[ind],PER_Y[ind],t='l')


#frequency
t_t0=t-min(t)
t_t0[1:5]/dt

dt=dt_select
lam_l=2*pi*(1:N)/(dt*N)
lam_l=matrix(lam_l,ncol=1,nrow=length(lam_l))
ind=1:(N/2)


PER_res_AM_FM_2018=apply(lam_l,1,periodogram,x=res_AM_FM_2018,times=t_t0)
plot(lam_l[ind],PER_res_AM_FM_2018[ind],t='l')

PER_Y=apply(lam_l,1,periodogram,x=Y,times=t_t0)
plot(lam_l[ind],PER_Y[ind],t='l')

W_y=apply(lam_l,1,W_fun,(t-min(t)))
plot(lam_l,W_y,t='l')

aux1=DFT_dar(W_y)$DFT_x
aux2=DFT_dar(PER_res_AM_FM_2018)$DFT_x
aux3=DFT_dar(PER_Y)$DFT_x

aux=Re((aux2*N_i)/(aux1*2*pi))
P_res_AM_FM_2018=Re(IN_DFT_dar(aux)$DFT_x)

aux=Re((aux3*N_i)/(aux1*2*pi))
P_Y=Re(IN_DFT_dar(aux)$DFT_x)

#####
plot(lam_l[1:1100],Mod(fft(Y))[1:1100]^2,t='h')
plot(lam_l[1:1100],Mod(fft(res_AM_FM_2018))[1:1100]^2,t='h')

####

pdf("plot1.pdf", width=10, height=7)

layout(matrix(1:6, 3,2, byrow = TRUE),widths=c(1,1), heights=c(1,1))

par(mar = c(3, 3.7, 0.5, 0.5))
plot(lam_l[ind],PER_Y[ind],t='h',ylab='',xlab='Frequency')
title(ylab='Periodogram of Y', line=2, font=2)
title(sub='Frequency', adj=0.5, line=1.9, font=2)

par(mar = c(3, 3.7, 0.5, 0.5))
plot(lam_l[ind],PER_res_AM_FM_2018[ind],t='h',ylab='',xlab='Frequency')
title(ylab='Periodogram of the residuals', line=2, font=2)
title(sub='Frequency', adj=0.5, line=1.9, font=2)

plot(lam_l[1:1100],Mod(fft(Y)^2)[1:1100],t='h',ylab='',xlab='')

plot(lam_l[10:1100],Mod(fft(res_AM_FM_2018)^2)[10:1100],t='h')

spectrum(Y,method=c("pgram"),log='no')

spectrum(res_AM_FM_2018,log='no')

dev.off()

Mod(fft(Y)^2)[1:5]
PER_Y[1:5]


pdf("plot2.pdf", width=10, height=7)

layout(matrix(1:4, 2,2, byrow = FALSE),widths=c(1,1), heights=c(1,1))

par(mar = c(3, 3.7, 0.5, 0.5))
plot(lam_l[ind],P_Y[ind],t='h',ylab='',xlab='Frequency')
title(ylab='PSD of Y', line=2, font=2)
title(sub='Frequency', adj=0.5, line=1.9, font=2)

par(mar = c(3, 3.7, 0.5, 0.5))
plot(lam_l[ind],P_res_AM_FM_2018[ind],t='h',ylab='',xlab='Frequency')
title(ylab='PSD of the residuals', line=2, font=2)
title(sub='Frequency', adj=0.5, line=1.9, font=2)

plot(lam_l[1:1100],Mod(fft(Y))[1:1100]^2,t='h',ylab='',xlab='')

plot(lam_l[10:1100],Mod(fft(res_AM_FM_2018))[10:1100]^2,t='h')


dev.off()


pdf("plot3.pdf", width=10, height=7)

par(mar = c(3, 3.7, 0.5, 0.5))
plot(lam_l[ind],W_y[ind],t='h',ylab='',xlab='Frequency')
title(ylab='Spectral window', line=2, font=2)
title(sub='Frequency', adj=0.5, line=1.9, font=2)

dev.off()



pdf("plot1.pdf", width=10, height=7)

layout(matrix(1:2, 2,1, byrow = FALSE),widths=c(1,1), heights=c(1,1))

plot(lam_l[1:1100],Mod(fft(Y))[1:1100]^2,t='h',ylab='',xlab='')

plot(lam_l[10:1100],Mod(fft(res_AM_FM_2018))[10:1100]^2,t='h')










#########################################################################################################


###MODEL 2 with B-SPLINES
setwd("~/Doctorado/Tesis/Modulation models/Paper ApJ 2021/R code Apj 2021 WN errors/Run profe Motta SAMPLE 2")
AIC_tot<-read.delim('AIC_V783_sample2.txt', header =FALSE, sep = " ")
dim(AIC_tot)#3840
head(AIC_tot)

ind_40=which(AIC_tot[,1]==40)

AIC_tot=AIC_tot[-ind_40,]

which.min(AIC_tot[,7])#2387
AIC_tot[2387,]
#30  0 0.1 0.1  0 10 2.5114220733e-06 0.002907769798


setwd("~/Doctorado/Tesis/Modulation models/Paper ApJ 2021/R code Apj 2021 WN errors")

K=19
K_benko=15
#new_freqs=c(18.3258794594958,19.9364443996956,21.5477609477191,23.1587514192813)

fe=0.6035
new_freqs=(11:14)*f0+fe
round(new_freqs,4)
f=c((1:K_benko)*f0,new_freqs)
w_hat=2*pi*f
p=1
bdeg=3
i=2387
lam_aux=c(AIC_tot[i,2],rep(AIC_tot[i,3],10),rep(AIC_tot[i,4],10),
          rep(AIC_tot[i,5],10),rep(AIC_tot[i,6],8))
nodos=AIC_tot[i,1]
fit=psfit3(x=t,xl=min(t)-0.00001,xr=max(t)+0.00001,y=Y,w1=w_hat,K=K,pord=p,ndx=nodos,bdeg=bdeg,lam=lam_aux)
ncol(fit$B)*(2*K+1)


hat_m_g=fit$hat_m_g
Y_hat_Bspline=fit$f.hat
res_Bspline=Y-Y_hat_Bspline
coeff=fit$theta
S_lam=fit$S_lambda
sigma2=fit$sigma2hat
hat_m_g=fit$hat_m_g
sen_m_g=fit$sd_m_g
error_train_Bspline=mean((Y-Y_hat_Bspline)^2)
error_train_Bspline

n=19
ind=seq(2,(2*n),by=2)
m_hat=hat_m_g[,1]
g1_hat=hat_m_g[,ind]
g2_hat=hat_m_g[,-c(1,ind)]  

m_sen=sen_m_g[,1]
g1_sen=sen_m_g[,ind]
g2_sen=sen_m_g[,-c(1,ind)]   


plot(t,res_Bspline,t='l')


layout(matrix(1:2, 2, 1, byrow = TRUE),widths=c(1,1), heights=c(1,1))
par(mar = c(3, 3.7, 0.6, 0))
plot(t,Y,ylim=c(max(Y),min(Y)),t='l',lwd=2)
lines(t,Y_hat_AM_FM_2018, col=2, lwd=1)
par(mar = c(3, 3.7, 0.6, 0))
plot(t,Y,ylim=c(max(Y),min(Y)),t='l',lwd=2)
lines(t,Y_hat_Bspline, col=2, lwd=1)


#ERRORS
round(c(error_train_AM_FM_2018),6)
#8e-06
round(c(error_train_Bspline),6)
#1e-06

library(lomb)
per=lsp(res_Bspline,times=t,type='frequency',ofac=1)
per$peak.at

per=lsp(res_AM_FM_2018,times=t,type='frequency',ofac=1)
per$peak.at

per=lsp(res_Bspline,times=t,from=0,to=1/(2*dt),type='frequency',ofac=1)
PER_Bspline=per$power
FREQ_Bspline=per$scanned
SIG_Bspline=per$sig.level

per=lsp(res_AM_FM_2018,times=t,from=0,to=2/(2*dt),type='frequency',ofac=1)
per=lsp(Y,times=t,from=0,to=2/(2*dt),type='frequency',ofac=1)
per=lsp(yy,times=tt,from=0,to=2/(2*dt),type='frequency',ofac=1)

#pdf("fit_benko.pdf", width=10, height=5)
par(mar = c(3, 3.7, 0.5, 0.5))
plot(per$scanned,per$power,t='h',ylab='',xlab='Frequency')
abline(h=per$sig.level,col=4,lty=2)
title(ylab='LS Periodogram of residuals', line=2, font=2)
title(sub='Frequency', adj=0.5, line=1.9, font=2)
#dev.off()

abline(v=47.325756,col=2)
f_star=47.325756
2*F_n-f_star
f_0

PER_BENKO=per$power
FREQ_BENKO=per$scanned
SIG_BENKO=per$sig.level

per=lsp(Y,times=t,from=0,to=1/(2*dt),type='frequency',ofac=1)
PER_DATA=per$power
FREQ_DATA=per$scanned
SIG_DATA=per$sig.level

par(mar = c(3, 3.5, 1, 1))
plot(FREQ_DATA,PER_DATA,t='h',xlab='',ylab='',lwd=2)
for(j in 1:15){
  abline(v=j*f0,col=2)
}
fe=0.6048
for(j in 1:14){
  abline(v=j*f0+fe,col=3)
}





pdf("V783_LS_PER_1.pdf", width=15, height=11)

maxx=max(PER_Bspline,PER_BENKO,PER_DATA)
minn=min(PER_Bspline,PER_BENKO,PER_DATA)

layout(matrix(1:3, 3,1, byrow = FALSE),widths=c(1,1), heights=c(1,1))

par(mar = c(3, 3.5, 1, 1))
plot(FREQ_DATA,PER_DATA,t='h',ylim=c(minn,maxx),xlab='',ylab='')
text(max(FREQ_DATA)/2,maxx-0.0003,cex=1.3,expression("LS periodogram of V783"))
title(sub=bquote(paste('Frequency [',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote('LS periodogram'), line=2, font=2)
for(j in 1:15){
  abline(v=j*f0,col=2)
}
fe=0.6048
for(j in 1:14){
  abline(v=j*f0+fe,col=3)
}


par(mar = c(3, 3.5, 1, 1))
plot(FREQ_Bspline,PER_Bspline,t='h',ylim=c(minn,maxx),xlab='',ylab='')
text(max(FREQ_Bspline)/2,maxx-0.0003,cex=1.3,expression("LS periodogram of the residuals obtained with time-varying parameters"))
title(sub=bquote(paste('Frequency [',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote('LS periodogram'), line=2, font=2)
for(j in 1:15){
  abline(v=j*f0,col=2)
}
fe=0.6048
for(j in 1:14){
  abline(v=j*f0+fe,col=3)
}

par(mar = c(3, 3.5, 1, 1))
plot(FREQ_BENKO,PER_BENKO,t='h',ylim=c(minn,maxx),xlab='',ylab='')
abline(v=23.1588,col=2)
abline(v=21.5478,col=2)
abline(v=19.9364,col=2)
abline(v=18.3259,col=2)
text(max(FREQ_BENKO)/2,maxx-0.0003,cex=1.3,expression("LS periodogram of the residuals obtained with time-invariant parameters"))
title(sub=bquote(paste('Frequency [',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote('LS periodogram'), line=2, font=2)
for(j in 1:15){
  abline(v=j*f0,col=2)
}
fe=0.6048
for(j in 1:14){
  abline(v=j*f0+fe,col=3)
}


dev.off()



############# HAT PSD of residuals
maxx=max(res_AM_FM_2018,res_Bspline)
minn=min(res_AM_FM_2018,res_Bspline)
plot(t,res_AM_FM_2018,t='l',ylim=c(minn,maxx))
plot(t,res_Bspline,t='l',ylim=c(minn,maxx))

#frequency
t_t0=t-min(t)
dt=dt_select
lam_l=2*pi*(1:N)/(dt*N)
lam_l=matrix(lam_l,ncol=1,nrow=length(lam_l))

PER_res_AM_FM_2018=apply(lam_l,1,periodogram,x=res_AM_FM_2018,times=t_t0)
plot(lam_l,PER_res_AM_FM_2018,t='l')

PER_res_Bspline=apply(lam_l,1,periodogram,x=res_Bspline,times=t_t0)
plot(lam_l,PER_res_Bspline,t='l')


s_PER_res_AM_FM_2018=ksmooth(lam_l, PER_res_AM_FM_2018, kernel ="normal", bandwidth = 7.2)$y
plot(lam_l,PER_res_AM_FM_2018,t='l')
lines(lam_l,s_PER_res_AM_FM_2018,col=2)

s_PER_res_Bspline=ksmooth(lam_l, PER_res_Bspline, kernel ="normal", bandwidth = 7.2)$y
plot(lam_l,PER_res_Bspline,t='l')
lines(lam_l,s_PER_res_Bspline,col=2)

plot(lam_l,s_PER_res_AM_FM_2018,t='l')
plot(lam_l,s_PER_res_Bspline,t='l')


W_y=apply(lam_l,1,W_fun,(t-min(t)))
plot(lam_l,W_y,t='l')

aux1=DFT_dar(W_y)$DFT_x
aux2=DFT_dar(s_PER_res_AM_FM_2018)$DFT_x

plot(lam_l,Re(aux1),t='l')
plot(lam_l,Re(aux2),t='l')
plot(lam_l,Re(aux2)/Re(aux1),t='l')

aux=Re((aux2*N_i)/(aux1*2*pi))
P_res_AM_FM_2018=Re(IN_DFT_dar(aux)$DFT_x)

aux=Re((DFT_dar(s_PER_res_Bspline)$DFT_x*N)/(DFT_dar(W_y)$DFT_x*2*pi))
P_res_Bspline=Re(IN_DFT_dar(aux)$DFT_x)

maxx=max(P_res_AM_FM_2018,P_res_Bspline)
minn=min(P_res_AM_FM_2018,P_res_Bspline)

plot(lam_l,P_res_AM_FM_2018,t='l',ylim=c(minn,maxx))
plot(lam_l,P_res_Bspline,t='l',ylim=c(minn,maxx))

sum(P_res_AM_FM_2018>0)
sum(P_res_Bspline>0)


ind=1:(N/2)

pdf("V783_fit_AIC_SAMPLE2.pdf", width=15, height=11)

layout(matrix(1:8, 4,2, byrow = FALSE),widths=c(1,1), heights=c(1,1))
#T-V MODEL

par(mar = c(3, 3.5, 1, 1))
minn_y=min(Y)-0.2
maxx_y=max(Y)+0.1
plot(t,Y,ylim=c(maxx_y,minn_y),t='l',lwd=2,col=2,ylab='',xlab='',xaxt='none')
lines(t,Y_hat_Bspline, col=1, lwd=1)
rug(t,col='grey')
title(sub='Time [BJD-2454833 d]', adj=0.5, line=0.5, font=2)
title(ylab='Brightness [mag]', line=2, font=2)
text(min(t)+(max(t)-min(t))/2,minn_y+0.03,cex=1.3,bquote('Fit obtained with time-varying parameters (our novel model)'))
legend("topright", legend=c('Observations','Fit'),col=c(2,"black"), lty=c(1,1), cex=1)

par(mar = c(3, 3.5, 1, 1))
minn_res=min(res_AM_FM_2018,res_Bspline)
maxx_res=max(res_AM_FM_2018,res_Bspline)
plot(t,res_Bspline,ylim=c(maxx_res,minn_res),t='l',ylab='',xlab='')
rug(t,col='grey')
title(sub='Time [BJD-2454833 d]', adj=0.5, line=2, font=2)
title(ylab='Resdiduals [mag]', line=2, font=2)
text(min(t)+(max(t)-min(t))/2,minn_res+0.001,cex=1.3,bquote('Residual of fit obtained with time-varying parameters (our novel model)'))

maxx=max(log10(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))+0.2
minn=min(log10(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))

par(mar = c(3, 3.5, 1, 1))
plot(log10(lam_l[ind]/(2*pi)),log10(P_res_Bspline[ind]),t='l',ylim=c(minn,maxx),
     xlab='',ylab='')
text(-0.15,maxx-0.2,cex=1.3,expression("log10 of estimated PSD of the residuals obtained with time-varying parameters (our novel model)"))
title(sub=bquote(paste('log10 of frequency [log10 ',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote(paste('log10 of ',hat(PSD))), line=2, font=2)

maxx=max(sqrt(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))+0.0008
minn=min(sqrt(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))

par(mar = c(3, 3.5, 1, 1))
plot(lam_l[ind]/(2*pi),sqrt(P_res_Bspline[ind]),t='l',ylim=c(minn,maxx),xlab='',ylab='')
text(max(lam_l[ind]/(2*pi))/2,maxx-0.0003,cex=1.3,expression("Square root of estimated PSD of the residuals obtained with time-varying parameters (our novel model)"))
title(sub=bquote(paste('Frequency [',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote(sqrt(hat(PSD))), line=2, font=2)


#BENKO MODEL

par(mar = c(3, 3.5, 1, 1))
plot(t,Y,ylim=c(maxx_y,minn_y),t='l',lwd=2,col=2,ylab='',xlab='',xaxt='none',yaxt='none')
lines(t,Y_hat_AM_FM_2018, col=1, lwd=1)
rug(t,col='grey')
title(sub='Time [BJD-2454833 d]', adj=0.5, line=0.5, font=2)
title(ylab='Brightness [mag]', line=1, font=2)
text(min(t)+(max(t)-min(t))/2,minn_y+0.03,cex=1.3,expression(paste('Fit obtained with time-invariant parameters (Benk',"\u00f6",' 2018)')))
#bquote('Fit obtained with Model 1'))
legend("topright", legend=c('Observations','Fit'),col=c(2,"black"), lty=c(1,1), cex=1)

par(mar = c(3, 3.5, 1, 1))
plot(t,res_AM_FM_2018,ylim=c(maxx_res,minn_res),t='l',ylab='',xlab='',yaxt='none')
rug(t,col='grey')
title(sub='Time [BJD-2454833 d]', adj=0.5, line=2, font=2)
title(ylab='Resdiduals [mag]', line=1, font=2)
text(min(t)+(max(t)-min(t))/2,minn_res+0.001,cex=1.3,bquote('Residual of fit obtained with time-invariant parameters'))

maxx=max(log10(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))+0.2
minn=min(log10(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))

par(mar = c(3, 3.5, 1, 1))
plot(log10(lam_l[ind]/(2*pi)),log10(P_res_AM_FM_2018[ind]),t='l',ylim=c(minn,maxx),
     xlab='',ylab='',yaxt='none')
text(-0.15,maxx-0.2,cex=1.3,expression("log10 of estimated PSD of the residuals obtained with time-invariant parameters"))
title(sub=bquote(paste('log10 of frequency [log10 ',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote(paste('log10 of ',hat(PSD))), line=1, font=2)


maxx=max(sqrt(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))+0.0008
minn=min(sqrt(c(P_res_AM_FM_2018[ind],P_res_Bspline[ind])))

par(mar = c(3, 3.5, 1, 1))
plot(lam_l[ind]/(2*pi),sqrt(P_res_AM_FM_2018[ind]),t='l',ylim=c(minn,maxx),xlab='',ylab='',yaxt='none')
text(max(lam_l[ind]/(2*pi))/2,maxx-0.0003,cex=1.3,expression("Square root of estimated PSD of the residuals obtained with time-invariant parameters"))
title(sub=bquote(paste('Frequency [',d^{-1},']' )), adj=0.5, line=2.0, font=2)
title(ylab=bquote(sqrt(hat(PSD))), line=1, font=2)
#axis(1, at=new_freqs[1],labels=bquote(paste("f'"[11])), col.axis="red", las=1)
#axis(1, at=new_freqs[2],labels=bquote(paste("f'"[12])), col.axis="red", las=1)
#axis(1, at=new_freqs[3],labels=bquote(paste("f'"[13])), col.axis="red", las=1)
#axis(1, at=new_freqs[4],labels=bquote(paste("f'"[14])), col.axis="red", las=1,srt = 10)
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

plot(0,0)
text(x=0.25,y=-0.9,labels = bquote(paste("f"[11])))
text(x=0.25,y=-0.85,labels = bquote(paste("'")))
axis(1, at=0.25,labels = FALSE, las=1,pos=-1)


pdf("V783_fit_AIC_g.pdf", width=14, height=10)

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
title(ylab=bquote(paste(hat(m),' and ',hat(u))), line=(nline+1.5), font=2)


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
  title(ylab=bquote(paste(hat(g)[1*','*~.(assay)],' and ',hat(h)[1*','*~.(assay)])), line=(nline+1.5), font=2)
  
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
  title(ylab=bquote(paste(hat(g)[1*','*~.(assay)],' and ',hat(h)[1*','*~.(assay)])), line=nline, font=2)
  

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
  title(ylab=bquote(paste(hat(g)[2*','*~.(assay)],' and ',hat(h)[2*','*~.(assay)])), line=nline, font=2)
  
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
  title(ylab=bquote(paste(hat(g)[2*','*~.(assay)],' and ',hat(h)[2*','*~.(assay)])), line=nline, font=2)
  
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
  title(ylab=bquote(paste(hat(g)[2*','*~.(assay)],' and ',hat(h)[2*','*~.(assay)])), line=nline, font=2)

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
title(ylab=bquote(paste(hat(g)," '"[1*','*~.(assay)])), line=(nline+1.5), font=2)


par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g1_hat[index,i+1]+1.96*g1_sen[index,i+1],g1_hat[index,i+1]-1.96*g1_sen[index,i+1]))
minn=min(cbind(g1_hat[index,i+1]+1.96*g1_sen[index,i+1],g1_hat[index,i+1]-1.96*g1_sen[index,i+1]))
plot(t,g1_hat[,i+1],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g1_hat[,i+1]+1.96*g1_sen[,i+1],lty=2)
lines(t,g1_hat[,i+1]-1.96*g1_sen[,i+1],lty=2)
rug(t,col='grey')
assay <- as.character(i+1)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[1*','*~.(assay)])), line=(nline), font=2)

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g1_hat[index,i+2]+1.96*g1_sen[index,i+2],g1_hat[index,i+2]-1.96*g1_sen[index,i+2]))
minn=min(cbind(g1_hat[index,i+2]+1.96*g1_sen[index,i+2],g1_hat[index,i+2]-1.96*g1_sen[index,i+2]))
plot(t,g1_hat[,i+2],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g1_hat[,i+2]+1.96*g1_sen[,i+2],lty=2)
lines(t,g1_hat[,i+2]-1.96*g1_sen[,i+2],lty=2)
rug(t,col='grey')
assay <- as.character(i+2)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[1*','*~.(assay)])), line=(nline), font=2)

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g2_hat[index,i]+1.96*g2_sen[index,i],g2_hat[index,i]-1.96*g2_sen[index,i]))
minn=min(cbind(g2_hat[index,i]+1.96*g2_sen[index,i],g2_hat[index,i]-1.96*g2_sen[index,i]))
plot(t,g2_hat[,i],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g2_hat[,i]+1.96*g2_sen[,i],lty=2)
lines(t,g2_hat[,i]-1.96*g2_sen[,i],lty=2)
rug(t,col='grey')
assay <- as.character(i)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[2*','*~.(assay)])), line=(nline), font=2)

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g2_hat[index,i+1]+1.96*g2_sen[index,i+1],g2_hat[index,i+1]-1.96*g2_sen[index,i+1]))
minn=min(cbind(g2_hat[index,i+1]+1.96*g2_sen[index,i+1],g2_hat[index,i+1]-1.96*g2_sen[index,i+1]))
plot(t,g2_hat[,i+1],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g2_hat[,i+1]+1.96*g2_sen[,i+1],lty=2)
lines(t,g2_hat[,i+1]-1.96*g2_sen[,i+1],lty=2)
rug(t,col='grey')
assay <- as.character(i+1)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[2*','*~.(assay)])), line=(nline), font=2)

par(mar = c(3, 3.7, 0.6, 0))
maxx=max(cbind(g2_hat[index,i+2]+1.96*g2_sen[index,i+2],g2_hat[index,i+2]-1.96*g2_sen[index,i+2]))
minn=min(cbind(g2_hat[index,i+2]+1.96*g2_sen[index,i+2],g2_hat[index,i+2]-1.96*g2_sen[index,i+2]))
plot(t,g2_hat[,i+2],t='l',ylab='',xlab='',ylim=c(maxx,minn),yaxt='none',lwd=1.5)
lines(t,g2_hat[,i+2]+1.96*g2_sen[,i+2],lty=2)
lines(t,g2_hat[,i+2]-1.96*g2_sen[,i+2],lty=2)
rug(t,col='grey')
assay <- as.character(i+2)
title(sub="BJD-2450000 [d]", adj=0.5, line=2, font=2)
title(ylab=bquote(paste(hat(g)," '"[2*','*~.(assay)])), line=(nline), font=2)

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
title(ylab=bquote(paste(hat(g)," '"[1*','*~.(assay)])), line=(nline+1.5), font=2)

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
title(ylab=bquote(paste(hat(g)," '"[2*','*~.(assay)])), line=(nline), font=2)

par(mar = c(3, 3.7, 0.6, 0))
plot.new()
par(mar = c(3, 3.7, 0.6, 0))
plot.new()

dev.off()


max(t)-min(t)
