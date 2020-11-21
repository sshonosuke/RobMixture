

# install packages 
library(MASS)
library(mvtnorm)
library(e1071)
library(mvnfast)
library(otrimle)
library(tclust)


###  RGM: Robust Gaussian Mixtrue  ###
# Y: (n, p)-matrix of observations 
# K: number of components
# gam: value of gamma in the weight
# C: eigen constraint
# maxitr: maximum number of iteration 
# alpha: tuning parameter outlier detection
# init: algorithm to compute initial values
RGM <- function(Y, K=2, gam=0.2, C=10, maxitr=1000, alpha=0.01, init="TCL"){
  EPS <- 10^(-4)   # torelence rate
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  
  # function for eigen-constraint
  Trunc <- function(xx, mm){
    xx[xx<mm] <- mm
    xx[xx>(C*mm)] <- C*mm
    return(xx)
  }
  
  # initial values
  PP <- rep(1/K, K)
  if(init=="TCL"){  
    fit <- tclust(Y, k=K)
    Mu <- fit$centers
  }
  if(init=="IML"){  
    fit <- otrimle(Y, G=K, monitor=F) 
    Mu <- fit$mean
  }
  Sig <- fit$cov
  for(k in 1:K){ Sig[,,k] <- 0.5*(Sig[,,k]+t(Sig[,,k])) }
  
  ### EEE iteration
  uu <- NA
  val <- c()
  
  for(j in 1:maxitr){
    Mu0 <- Mu
    
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
    val[j] <- sqrt(sum((Mu-Mu0)^2))/sum(abs(Mu0))
    if( val[j] < EPS ){ break }
    if( j>1 & min(abs(val[j]-val[1:(j-1)])) < 10^(-7) ){ break }
  }
  
  # outlier detection
  group <- apply(uu, 1, which.max)
  out <- list()
  for(k in 1:K){
    ID <- (1:N)[group==k]
    Nk <- length(ID)
    if(Nk>0){
      MC <- 10000 
      dens <- dmvnorm(Y, Mu[,k], Sig[,,k])[group==k]
      ref.rn <- mvrnorm(MC, Mu[,k], Sig[,,k])
      ref.val <- dmvnorm(ref.rn, Mu[,k], Sig[,,k])
      th <- quantile(ref.val, prob=alpha)
      out[[k]] <- ID[dens<th]
    }
  }
  out <- unlist(out)
  if(gam>0){ group[out] <- 0 }
  
  # BIC
  Dens <- matrix(NA, N, K)
  for(k in 1:K){
    Dens[,k] <- dmvnorm(Y, Mu[,k], Sig[,,k])
  }
  dens <- apply(t(Dens)*PP, 2, sum)    # estimated density 
  ML <- N*mean( log(dens[-out]) )
  BIC <- -2*ML + log(N)*(p*K+K*p*(p+1)/2+K-1)
  
  # Summary
  Res <- list(Mu=Mu, Sigma=Sig, Pi=PP, out=out, BIC=BIC, cluster=group, post.prob=uu, itr=j)
  return(Res)
}










###  Robust selection of the number of components in RGM  ###
# Y: (n, p)-matrix of observations 
# K.min: minimum number of components
# K.max: maximum number of components
# gam: value of gamma in the weight
# C: eigen constraint
# maxitr: maximum number of iteration 
RGM.select <- function(Y, K.min=2, K.max=7, gam=0.2, C=10, maxitr=100, alpha=0.01){
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  K.set <- K.min:K.max
  L <- length(K.set)
  
  BIC <- c()
  Fit <- list()
  for(l in 1:L){
    # Fitting
    sK <- K.set[l]
    fit <- RGM(Y, K=sK, gam=gam, C=C, maxitr=maxitr, alpha=alpha, init="IML")
    BIC[l] <- fit$BIC
    Fit[[l]] <- fit
  }
  
  # optimal number of components
  opt.K <- K.set[which.min(BIC)]
  opt.fit <- Fit[[which.min(BIC)]]
  Res <- list(opt.K=opt.K, opt.fit=opt.fit, BIC=BIC)
  return(Res)
}




