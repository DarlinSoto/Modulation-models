##########################
# AUTHOR: DARLIN SOTO    #
# DATE: OCTOBER 15, 2021 #
# EMAIL: dmsoto1@uc.cl   #
##########################

#This code reproduces Figure 1 in the paper Periodic Variable Stars Modulated by Time-Varying Parameters.

library(splines)

n=500
set.seed(1234) # for reproducible results


#SIMULATION DATA
sigma2=1
x = sort(round(runif(n,0,55),2))
w=0.2*pi
trend=-0.05*x
g11=-0.0002*x+0.0003*x^2
g21=1-0.0005*x
fX = trend+g11*cos(w*x)+g21*sin(w*x)
z <-rnorm(n, mean = 0, sd = sqrt(sigma2))
Y = fX + z 


#PARAMETERS B-SPLINE
num_knots=10
pordd=1
tau=seq(0,200,by=5)
tau_length=length(tau)
tau=matrix(tau,ncol=1,nrow=tau_length)


###############################################

AIC=matrix(0,ncol=1,nrow=tau_length)
Y_hat=matrix(0,ncol=tau_length,nrow=n)

for(l in 1:length(tau)){
  
  fit=psfit(x=x,xl=min(x)-0.0000001,xr=max(x)+0.0000001,y=Y,w1=w,K=1,pord=pordd,
             ndx=num_knots,bdeg=3,lam=c(tau[l],tau[l],tau[l]))
  
  AIC[l]=fit$AIC
  Y_hat[,l]=fit$f.hat
}

Y_WN=Y
x_WN=x
fX_WN=fX
Y_hat_WN=Y_hat
AIC_WN=AIC
tau_WN=tau
z_WN=z
tau_length_WN=tau_length


pdf("PSE_WN.pdf",width=10,height=3)

layout(matrix(c(1,1,2), 1,3, byrow = TRUE),widths=c(1,1), heights=c(1,1))

maxx=max(Y_WN)
minn=min(Y_WN)-0.5
par(mar = c(3, 3.5, 1, 1))
plot(x_WN,Y_WN,col="grey",pch=16,ylab='',xlab='',ylim=c(minn,maxx))
lines(x_WN,Y_WN-z_WN,t='l',lwd=1.5)
lines(x_WN,Y_hat_WN[,1],col="orange",lty=2)
lines(x_WN,Y_hat_WN[,which.min(AIC_WN)],col="blue",lty=2)
lines(x_WN,Y_hat_WN[,length(tau_WN)],col="green",lty=2)
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(paste('y, ', mu, ' and ', hat(y))), line=2, font=2)
legend("topright", 
       legend=c((expression(paste(tau," = 0"))),(expression(paste(tau," = 30"))),(expression(paste(tau," = 200")))),
       col=c("orange","blue","green"), lty=2, cex=1)
rug(x_WN,col='grey')

maxx=max(AIC_WN)
minn=min(AIC_WN)
par(mar = c(3, 3.5, 1, 1))
plot(tau_WN,AIC_WN,t='l',col="black",lwd=1,ylab='',xlab='',ylim=c(minn,maxx))
points(tau_WN[c(1,which.min(AIC_WN),tau_length_WN)],AIC_WN[c(1,which.min(AIC_WN),tau_length_WN)],pch=19,col=c("orange","blue","green"))
title(sub=bquote(tau), adj=0.5, line=2, font=2)
title(ylab=bquote(AIC), line=2, font=2)

dev.off()





