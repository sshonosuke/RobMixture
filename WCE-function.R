###-----------------------------------------------------------###
###    Function for robust mixture modeling using WCE         ###
###-----------------------------------------------------------###

## This code includes the following 4 functinons: 
#  RGM: robust Gaussian mixture
#  RGM.selection: robust selection of the number of components in Gaussian mixture
#  RSNM: robust skew normal mixture
#  RMOE: robust mixture of experts (two component mixture of linear regressions)



# install packages 
library(MASS)
library(mvtnorm)
library(e1071)
library(mvnfast)
library(otrimle)



###  RGM: Robust Gaussian Mixtrue  ###
# Y: (n, p)-matrix of observations 
# K: number of components
# gam: value of gamma in the weight
# C: eigen constraint
# maxitr: maximum number of iteration 
RGM <- function(Y, K=2, gam=0.5, C=10, maxitr=1000){
  EPS <- 10^(-6)   # torelence rate
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  
  # function for eigen-constraint
  Trunc <- function(xx, mm){
    xx[xx<mm] <- mm
    xx[xx>(C*mm)] <- C*mm
    return(xx)
  }
  
  # initial values
  PP <- rep(1/K,K)
  fit <- otrimle(Y,G=K,erc=C,monitor=F)
  Mu <- fit$mean
  Sig <- fit$cov
  for(k in 1:K){ Sig[,,k] <- 0.5*(Sig[,,k]+t(Sig[,,k])) }
  
  ### EEE iteration
  val <- 0
  uu <- NA
  
  for(j in 1:maxitr){
    uu0 <- uu
    val0 <- val
    ##  E-step
    Dens <- matrix(NA,N,K)
    for ( k in 1:K ) {
      Dens[,k] <- dmvnorm(Y,Mu[,k],Sig[,,k])
    }
    uu <- t(t(Dens)*PP)
    uu <- uu/apply(uu,1,sum)
    if(K==1){ uu <- matrix(rep(1,N), N, 1) }
    ww <- Dens^(gam)
    
    B <- c()
    G <- c()
    for(k in 1:K){
      B[k] <- det(Sig[,,k])^(-gam/2)*(2*pi)^(-p*gam/2)*(1+gam)^(-p/2)
      G[k] <- gam*B[k]/(1+gam)
    }
    
    ##  EE-step
    PP <- apply(t(uu*ww)/B,1,sum) / sum(t(uu*ww)/B)   # mixing probability 
    D <- matrix(NA,p,K)
    H <- list()
    for(k in 1:K){
      # Mean
      Mu[,k] <- apply(uu[,k]*ww[,k]*Y,2,sum) / sum(uu[,k]*ww[,k])
      # Variance 
      resid <- t(Y)-Mu[,k]
      SS <- resid%*%diag(uu[,k]*ww[,k])%*%t(resid) / ( sum(uu[,k]*ww[,k])-sum(uu[,k])*G[k] )
      SV <- svd(SS)
      D[,k] <- SV$d
      H[[k]] <- SV$u
    }
    # eigen-constraint
    opt.M <- function(mm){
      tD <- Trunc(D,mm)
      sum( dim(Y)[[1]]*PP*apply(log(tD)+D/tD,2,sum) )
    }
    m.ast <- optim(par=1, fn=opt.M, method="L-BFGS-B", lower=0.001)$par
    
    # update Sigma
    for(k in 1:K){
      Sig[,,k] <- t(H[[k]])%*%diag(Trunc(D[,k],m.ast))%*%H[[k]]
      Sig[,,k] <- 0.5*(Sig[,,k]+t(Sig[,,k]))
    }
    
    # Convergence check
    Dens <- matrix(NA,N,K)
    for ( k in 1:K ) {
      Dens[,k] <- dmvnorm(Y, Mu[,k], Sig[,,k])
    }
    val <- sum( log(apply(t(Dens)*PP, 2, sum)) )
    if( is.na(val) ){ val=0 }
    if( abs(val)==Inf ){ val=1000 }
    if( abs(val-val0)/(abs(val0)+0.01) < EPS ){ break }
  }
  
  # Summary
  Res <- list(Mu=Mu, Sigma=Sig, Pi=PP)
  return(Res)
}






###  Robust selection of the number of components in RGM  ###
# Y: (n, p)-matrix of observations 
# K.min: minimum number of components
# K.max: maximum number of components
# gam: value of gamma in the weight
# C: eigen constraint
# maxitr: maximum number of iteration 
RGM.selection <- function(Y, K.min=2, K.max=7, gam=0.5, C=10, maxitr=1000){
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  K.set <- K.min:K.max
  L <- length(K.set)
  
  BIC <- c()
  for(l in 1:L){
    # Fitting
    sK <- K.set[l]
    fit <- RGM(Y, K=sK, gam=gam, C=C, maxitr=maxitr)
    
    #  BIC 
    Mu <- fit$Mu
    Sig <- fit$Sig
    PP <- fit$Pi
    Dens <- matrix(NA,N,sK)
    for(k in 1:sK){
      Dens[,k] <- dmvnorm(Y, Mu[,k], Sig[,,k])
    }
    dens <- apply(t(Dens)*PP, 2, sum)  # estimated density
    out.max <- N*0.15    # the maximum number of outliers
    ww <- dens^(gam)
    Ind <- (1:N)[ww<median(ww)/10]
    if(length(Ind)>out.max){
      sub <- order(ww[Ind], decreasing=F)[1:out.max]
      Ind <- Ind[sub]
    }
    ML <- N*mean(log(dens[-Ind]))
    BIC[l] <- -2*ML+log(N)*(p*sK+sK*p*(p+1)/2+sK-1)
  }
  
  # optimal number of components
  opt.K <- K.set[which.min(BIC)]
  return(opt.K)
}






###  RSNM: Robust Skew Normal Mixtrue  ###
# Y: (n, p)-matrix of observations 
# K: number of components
# gam: value of gamma in the weight
# C: eigen constraint
# maxitr: maximum number of iteration 

RSNM <- function(Y, K=2, gam=0.5, C=5, maxitr=1000){
  EPS <- 10^(-6)   # torelence rate
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  
  # skew normal density function 
  Density <- function(Mu, Sig, Psi){
    Om <- array(NA,c(p, p, K))
    som <- array(NA,c(p, p, K))
    Alpha <- matrix(NA, p, K)
    for(k in 1:K){
      Om[,,k] <- Sig[,,k]+Psi[,k]%*%t(Psi[,k])
      som[,,k] <- diag( sqrt( diag(Om[,,k]) ) )
      denom <- as.vector( sqrt( 1-Psi[,k]%*%solve(Om[,,k])%*%Psi[,k] ) )
      Alpha[,k] <- as.vector( som[,,k]%*%solve(Om[,,k])%*%Psi[,k] )/denom
    }
    f <- matrix(NA, N, K)
    for(k in 1:K){
      f[,k] <- 2*dmvnorm(Y,Mu[,k],Om[,,k])*pnorm( Alpha[,k]%*%solve(som[,,k])%*%(t(Y)-Mu[,k]) )
    }
    return(f)
  }
  
  # Initial values
  PP <- rep(1/K,K)
  Sig <- array(NA,c(p, p, K))
  for(k in 1:K){ Sig[,,k] <- diag(p) }
  fit.init <- otrimle(Y, G=K, erc=C, monitor=F)
  Mu <- fit.init$mean
  Psi <- matrix(NA,p,K)
  for(k in 1:K){
    Psi[,k] <- apply(Y[fit.init$cluster==k,], 2, skewness)
  }
  
  
  ##  EEE Iteration
  val <- 0
  for(j in 1:maxitr){
    val0 <- val
    ## E-step 
    Tau <- matrix(NA, N, K)
    Delta <- matrix(NA, N, K)
    TT <- matrix(NA, N, K)
    M <- matrix(NA, N, K)
    U <- matrix(NA, N, K)
    for(k in 1:K){
      c1 <- as.vector( t(Psi[,k])%*%solve(Sig[,,k])%*%Psi[,k] )
      c2 <- as.vector( t(Psi[,k])%*%solve(Sig[,,k])%*%(t(Y)-Mu[,k]) )
      c3 <- diag( t(t(Y)-Mu[,k])%*%solve(Sig[,,k])%*%(t(Y)-Mu[,k]) )
      Tau[,k] <- 1/sqrt( c1+1 )
      Delta[,k] <- c2/(c1+1)
      TT[,k] <- ( gam*c1 + 1/Tau[,k]^2 )^(-1/2)
      M[,k] <- TT[,k]^2*( gam*c2 + Delta[,k]/Tau[,k]^2)
      ratio <- M[,k]/TT[,k]; ratio[ratio<(-20)]=-20
      a1 <- det(Sig[,,k])^(-gam/2)*(2*pi)^(-gam*p/2)*pnorm(ratio)/pnorm(ratio)*(TT[,k]/Tau[,k])
      a2 <- exp(-0.5*gam*c3 - 0.5*Delta[,k]^2/Tau[,k]^2 + 0.5*ratio^2 )
      U[,k] <- a1*a2
    }
    Dens=Density(Mu,Sig,Psi)   # density values
    uu=t(t(Dens)*PP)
    uu=uu/apply(uu,1,sum)   # posteiror probability
    
    # calculation of V
    V0 <- U
    V1 <- matrix(NA,N,K)
    V2 <- matrix(NA,N,K)
    for(k in 1:K){
      ratio <- M[,k]/TT[,k]; ratio[ratio<(-20)]=-20
      V1[,k] <- U[,k]*( M[,k]+TT[,k]*dnorm(ratio)/pnorm(ratio) )
      V2[,k] <- U[,k]*( M[,k]^2+TT[,k]^2+M[,k]*TT[,k]*dnorm(ratio)/pnorm(ratio) )
    }
    
    # calculation of adjustment constants
    B <- c()
    GG <- c()
    for(k in 1:K){
      B[k] <- det(Sig[,,k])^(-gam/2)*(2*pi)^(-p*gam/2)*(1+gam)^(-p/2)
      GG[k] <- gam*B[k]/(1+gam)
    }
    
    ##  M-step
    PP <- apply(t(uu*V0)/B,1,sum)/sum(t(uu*V0)/B)
    # Mu
    for(k in 1:K){
      num <- apply( uu[,k]*(V0[,k]*Y-V1[,k]%*%t(Psi[,k])), 2,sum)
      denom <- sum(uu[,k]*V0[,k])
      Mu[,k] <- num/denom
    }
    # Psi
    for(k in 1:K){
      num <- apply( uu[,k]*V1[,k]*t(t(Y)-Mu[,k]), 2,sum)
      denom <- sum(uu[,k]*V2[,k])
      Psi[,k] <- num/denom
    }
    # Sigma (decomposition matrix)
    A <- array(NA,c(N, K, p, p))
    for(k in 1:K){
      resid <- t(t(Y)-Mu[,k])
      for(i in 1:N){
        A[i,k,,] <- V0[i,k]*resid[i,]%*%t(resid[i,])-2*V1[i,k]*resid[i,]%*%t(Psi[,k])+V2[i,k]*Psi[,k]%*%t(Psi[,k])
      }
    }
    
    D <- matrix(NA, p, K)
    H <- list()
    for(k in 1:K){
      num <- matrix(0, p, p)
      for(i in 1:N){
        num <- num+uu[i,k]*A[i,k,,]
      }
      denom <- sum(uu[,k]*V0[,k])-GG[k]*sum(uu[,k])
      SS <- num/denom
      SV <- svd(SS)
      D[,k] <- SV$d
      H[[k]] <- SV$u
    }
    
    # eigen-constraint
    Trunc <- function(xx,mm){
      xx[xx<mm] <- mm
      xx[xx>(C*mm)] <- C*mm
      return(xx)
    }
    opt.M <- function(mm){
      tD <- Trunc(D,mm)
      sum( dim(Y)[[1]]*PP*apply(log(tD)+D/tD,2,sum) )
    }
    m.ast <- optim(par=1,fn=opt.M,method="L-BFGS-B",lower=0.001)$par
    
    # update Sigma
    for(k in 1:K){
      Sig[,,k] <- t(H[[k]])%*%diag(Trunc(D[,k],m.ast))%*%H[[k]]
      Sig[,,k] <- 0.5*(Sig[,,k]+t(Sig[,,k]))
    }
    
    ## convergence check 
    Dens <- Density(Mu,Sig,Psi)
    val <- mean( log( apply((t(Dens)*PP),2,sum) ))
    if(is.na(val)){ val <- 0 }
    if(abs(val)==Inf){ val <- 1000 }
    if( abs(val-val0)/(abs(val0)+0.01) < EPS ){ break }
  }
  
  # Summary
  Res <- list(Mu=Mu, Sigma=Sig, Psi=Psi, Pi=PP)
  return(Res)
}






###  Robust Mixture-of-experts  ###
##  (Two component mixture of linear regression with a single covariate)  ##

# logistic function (used in the main function)
logistic <- function(x){ 
  val <- 1/(1+exp(-x)) 
  val[val<10^(-8)] <- 10^(-8)
  val[val>0.99999] <- 0.99999
  return(val)
}

##  main function  ##
RMOE <- function(Y, X, gam=0.2, maxitr=1000){
  K <- 2
  EPS <- 10^(-6)   # torelence rate
  N <- length(Y)
  
  # initial values
  PP <- rep(0.5, N)
  Beta <- matrix(c(2, 0, 0, 1), 2, 2)
  Eta <- c(0, 0)
  Sig <- rep(sd(Y), 2)
  XX <- cbind(1, X)
  
  ### EEE iteration
  for(j in 1:maxitr){
    Beta0 <- Beta
    # E-step
    PP <- logistic(Eta[1] + Eta[2]*X)
    log.dens1 <- dnorm(Y, Beta[1,1]+Beta[2,1]*X, Sig[1], log=T)
    log.dens2 <- dnorm(Y, Beta[1,2]+Beta[2,2]*X, Sig[2], log=T)
    log.uu <- cbind(log(PP)+log.dens1, log(1-PP)+log.dens2) 
    mm <- apply(log.uu, 1, max)
    uu <- exp(log.uu-mm)
    uu <- uu/apply(uu,1,sum)
    ww <- cbind(exp(log.dens1)^(gam), exp(log.dens2)^(gam))
    
    G <- c()
    for(k in 1:K){
      G[k] <- gam/(1+gam) * (Sig[k]^2)^(-gam/2)*(2*pi)^(-gam/2)*(1+gam)^(-1/2)
    }
    
    ## EE-step
    # regression coefficients and variance
    for(k in 1:K){
      if(sum(uu[,k])>3){
        Beta[,k] <- ginv(t(XX)%*%(uu[,k]*ww[,k]*XX))%*%t(XX)%*%(uu[,k]*ww[,k]*Y)
        mu <- Beta[1,k]+Beta[2,k]*X
        SS <- sum(uu[,k]*ww[,k]*(Y-mu)^2) / ( sum(uu[,k]*ww[,k])-sum(uu[,k])*G[k] )
        if(SS<0){
          SS <- sum(uu[,k]*ww[,k]*(Y-mu)^2) / sum(uu[,k]*ww[,k])
        }
        Sig[k] <- sqrt(SS)
      }
    }
    
    # eta in mixture proportion 
    log.dens1 <- dnorm(Y, Beta[1,1]+Beta[2,1]*X, Sig[1], log=T)
    log.dens2 <- dnorm(Y, Beta[1,2]+Beta[2,2]*X, Sig[2], log=T)
    ww <- cbind(exp(log.dens1)^(gam), exp(log.dens2)^(gam))
    B <- c()
    for(k in 1:K){
      B[k] <- (Sig[k]^2)^(-gam/2)*(2*pi)^(-gam/2)*(1+gam)^(-1/2)
    }
    opt <- function(eta){
      pp <- logistic(eta[1]+eta[2]*X)
      val <- sum( t( t(uu*ww)/B )* cbind(log(pp), log(1-pp)) )
      val <- max(val, -10^10)
      return(-val)
    }
    Eta <- optim(par=c(0,0), lower=c(-10,-10), upper=c(10,10), fn=opt, method="L-BFGS-B")$par
    
    # Convergence check
    dd <- sum(abs(Beta-Beta0)) / (sum(abs(Beta0))+0.001)
    if( dd < EPS ){ break }
  }
  
  # Summary
  Para <- list(Beta=Beta, Sigma=Sig, Eta=Eta)
  return(Para)
}
















