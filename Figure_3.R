##########################
# AUTHOR: DARLIN SOTO    #
# DATE: OCTOBER 16, 2021 #
# EMAIL: dmsoto1@uc.cl   #
##########################

library(NISTunits)
library(splines)
library(gtools)
library(Matrix)


#READ SIMULATED BLAZHKO STAR
Blazhko_sim=read.table("Blazhko_sim.txt",header=TRUE, sep = ",")
Y_t=Blazhko_sim[,1]
t=Blazhko_sim[,2]


#FIT
knots_selec=15
pordd=1
lam_m=5
lam_g11=1
lam_g12=0.1
lam_g13=0.1
lam_g14=0.1
lam_g21=0.1
lam_g22=0.1
lam_g23=1
lam_g24=4
lam_selec=c(lam_m,lam_g11,lam_g21,lam_g12,lam_g22,lam_g13,lam_g23,lam_g14,lam_g24)

fit=psfit(x=t,xl=min(t)-0.0001,xr=max(t)+0.0001,y=Y_t,w1=w,K=K,pord=pordd,ndx=knots_selec,bdeg=3,lam=lam_selec)

f_hat=fit$f.hat
resid=Y_t-f_hat
hat_m_g=fit$hat_m_g
sd=fit$sd_y
sigma2=fit$sigma2hat
sen_m_g=fit$sd_m_g

hat_m=hat_m_g[,1]
hat_g11=hat_m_g[,2]
hat_g21=hat_m_g[,3]
hat_g12=hat_m_g[,4]
hat_g22=hat_m_g[,5]
hat_g13=hat_m_g[,6]
hat_g23=hat_m_g[,7]
hat_g14=hat_m_g[,8]
hat_g24=hat_m_g[,9]

sen_m=sen_m_g[,1]
sen_g11=sen_m_g[,2]
sen_g21=sen_m_g[,3]
sen_g12=sen_m_g[,4]
sen_g22=sen_m_g[,5]
sen_g13=sen_m_g[,6]
sen_g23=sen_m_g[,7]
sen_g14=sen_m_g[,8]
sen_g24=sen_m_g[,9]



pdf("RRLyrae_sim_fit.pdf", width=10, height=11)

nline=1.8

layout(matrix(c(1,1,2,2,3,3,4:11), 7, 2, byrow = TRUE),widths=c(1,1), heights=c(1,1))

par(mar = c(3, 3.5, 0.5, 0))
plot(t,Y_t,col="grey",pch=16,ylab='',xlab='',ylim=c(max(Y_t),min(Y_t)),xaxt='none')
lines(t,mu,t='l',col=2,lwd=1.5)
lines(t,f_hat,t='l',col=1)
lines(t,f_hat+1.96*sd,t='l',lty=2,col=1)
lines(t,f_hat-1.96*sd,t='l',lty=2,col=1)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(paste('Y, ', mu, ' and ', hat(Y))), line=2, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,resid,t='l',ylim=c(min(resid),max(resid)),xlab='',ylab='',xaxt='none')
rug(t,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab='Residuals', line=2.1, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,mm,col=2,lwd=1.5,type='l',ylab='',xlab='',ylim=c(0.06,-0.03))
lines(t,hat_m,col=1,lwd=1)
lines(t,hat_m+1.96*sen_m,col=1,lwd=1,lty=2)
lines(t,hat_m-1.96*sen_m,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=1.5, font=2)
title(ylab=bquote(m*' and '* hat(m)), line=nline, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,g11,col=2,lwd=1.5,type='l',ylab='',xlab='',xaxt="none",ylim=c(0.12,-0.07))
lines(t,hat_g11,col=1,lwd=1)
lines(t,hat_g11+1.96*sen_g11,col=1,lwd=1,lty=2)
lines(t,hat_g11-1.96*sen_g11,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[1*','*1]*' and '* hat(g)[1*','*1]), line=nline, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,g21,col=2,lwd=1.5,type='l',ylab='',xlab='',xaxt="none",ylim=c(1,-0.2))
lines(t,hat_g21,col=1,lwd=1)
lines(t,hat_g21+1.96*sen_g21,col=1,lwd=1,lty=2)
lines(t,hat_g21-1.96*sen_g21,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[2*','*1]*' and '* hat(g)[2*','*1]), line=nline, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,g12,col=2,lwd=1.5,type='l',ylab='',xlab='',xaxt="none",ylim=c(0.26,-0.07))
lines(t,hat_g12,col=1,lwd=1)
lines(t,hat_g12+1.96*sen_g12,col=1,lwd=1,lty=2)
lines(t,hat_g12-1.96*sen_g12,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[1*','*2]*' and '* hat(g)[1*','*2]), line=nline, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,g22,col=2,lwd=1.5,type='l',ylab='',xlab='',xaxt="none",ylim=c(0.1,-0.35))
lines(t,hat_g22,col=1,lwd=1)
lines(t,hat_g22+1.96*sen_g22,col=1,lwd=1,lty=2)
lines(t,hat_g22-1.96*sen_g22,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[2*','*2]*' and '* hat(g)[2*','*2]), line=nline, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,g13,col=2,lwd=1.5,type='l',ylab='',xlab='',xaxt="none",ylim=c(0.1,-0.35))
lines(t,hat_g13,col=1,lwd=1)
lines(t,hat_g13+1.96*sen_g13,col=1,lwd=1,lty=2)
lines(t,hat_g13-1.96*sen_g13,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[1*','*3]*' and '* hat(g)[1*','*3]), line=nline, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,g23,col=2,lwd=1.5,type='l',ylab='',xlab='',xaxt="none",ylim=c(0.12,-0.07))
lines(t,hat_g23,col=1,lwd=1)
lines(t,hat_g23+1.96*sen_g23,col=1,lwd=1,lty=2)
lines(t,hat_g23-1.96*sen_g23,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=0.5, font=2)
title(ylab=bquote(g[2*','*3]*' and '* hat(g)[2*','*3]), line=nline, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,g14,col=2,lwd=1.5,type='l',ylab='',xlab='',ylim=c(0.25,-0.07))
lines(t,hat_g14,col=1,lwd=1)
lines(t,hat_g14+1.96*sen_g14,col=1,lwd=1,lty=2)
lines(t,hat_g14-1.96*sen_g14,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(g[1*','*4]*' and '* hat(g)[1*','*4]), line=nline, font=2)

par(mar = c(3, 3.5, 0.5, 0))
plot(t,g24,col=2,lwd=1.5,type='l',ylab='',xlab='',ylim=c(0.07,-0.03))
lines(t,hat_g24,col=1,lwd=1)
lines(t,hat_g24+1.96*sen_g24,col=1,lwd=1,lty=2)
lines(t,hat_g24-1.96*sen_g24,col=1,lwd=1,lty=2)
rug(t,col='grey')
title(sub="Time", adj=0.5, line=2, font=2)
title(ylab=bquote(g[2*','*4]*' and '* hat(g)[2*','*4]), line=nline, font=2)

dev.off()

