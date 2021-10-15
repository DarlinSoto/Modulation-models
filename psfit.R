##########################
# AUTHOR: DARLIN SOTO    #
# DATE: OCTOBER 15, 2021 #
# EMAIL: dmsoto1@uc.cl   #
##########################

psfit3 <- function(x,xl,xr,y,w1,K,pord,ndx,bdeg,lam){
  
  #B-spline
  n <- length(y)
  B <- bspline(x,xl,xr,ndx,bdeg=bdeg)$B
  knots<-bspline(x,xl,xr,ndx,bdeg=bdeg)$knots
  
  #Design matrix
  j=1:K
  C=list()
  S=list()
  CB=list()
  SB=list()
  X=cbind(B)
  for(k in 1:K){
    C[[k]]<-diag(cos(w1[k]*x))
    S[[k]]<-diag(sin(w1[k]*x))
    CB[[k]]=C[[k]]%*%B
    SB[[k]]=S[[k]]%*%B
    X=cbind(X,CB[[k]],SB[[k]])
  }
  nb <- ncol(B)
  nX <- ncol(X)
    
  # Construct penalty stuff
  D=diff(diag(nb), diff = pord)
  lambda=diag(lam)
  P <- kronecker(sqrt(lambda),D)
  
  # Fit
  f = lsfit(rbind(X, P), c(y, rep(0,nrow(P))), intercept = FALSE)
  h = hat(f$qr)[1:n]
  theta = f$coef
  
  ########################################
  
  #Estimation Y
  f.hat = X %*% theta
  
  #Estimation trend and amplitudes
  BB=B
  for(k in 1:K){
    aux<-cbind(B,B)
    BB=cbind(BB,aux)
  }
  Ncol_C=dim(B)[2]
  aux=rep(1,dim(B)[2])
  
  Diag_1=list()
  index=1:(2*K*Ncol_C+Ncol_C)
  for(k in 0:(2*K)){
    aux=index[(k*Ncol_C+1):((k+1)*Ncol_C)]
    Diag_1[[k+1]]=matrix(0,(2*K*Ncol_C+Ncol_C),(2*K*Ncol_C+Ncol_C))
    Diag_1[[k+1]][aux,aux]=diag(1,Ncol_C)
  }
  
  hat_m_g<-matrix(0,nrow=n,ncol=(2*K+1))
  
  for(k in 1:(2*K+1)){
    hat_m_g[,k]=BB%*%Diag_1[[k]]%*%theta
  }
  
  #Estimation sigma^2
  R_lambda=solve(t(X)%*%X+t(P)%*%P)%*%t(X)
  S=X%*%R_lambda
  
  rss = sum((y-f.hat)^2)
  df=sum(diag(S))
  df_res=sum(diag(S%*%t(S)))-2*df+n
  sigma2hat=rss/(n - df)
  sigma2hat2=rss/df_res
  
  ##################################
  
  #Variance of Y
  Covz2= sigma2hat * S%*%t(S)
  sd_y=sqrt(diag(Covz2))
  
  #Variance of trend and amplitudes
  sd_m_g=matrix(0,ncol=(2*K+1),nrow=n)
  for(k in 1:(2*K+1)){
    S_m_g=BB%*%Diag_1[[k]]%*%R_lambda
    Cov_m_g=sigma2hat * S_m_g%*%t(S_m_g)
    sd_m_g[,k]=sqrt(diag(Cov_m_g))
  }
  
  
  #Estimation with tau=0
  f = lsfit(X, y, intercept = FALSE)
  theta_0 = f$coef
  f.hat_0 = X %*% theta_0
  sigma2_0 <- mean((y-f.hat_0)^2)
  
  #AIC
  deviance=sum((y-f.hat)^2)
  AIC=deviance/n+2*df*sigma2_0/n
  
  output <- list(X=X,B=B,sigma2hat=sigma2hat,sigma2hat2=sigma2hat2,rss=rss,
                 theta=theta,f.hat=f.hat,hat_m_g=hat_m_g,
                 sd_y=sd_y,sd_m_g=sd_m_g,
                 AIC=AIC,BIC=BIC,R_lambda=R_lambda)
  return(output)
}