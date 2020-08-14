library(MASS)
library(mvtnorm)
library(e1071)
library(mvnfast)
library(otrimle)
library(truncnorm)


###  RSNM: Robust Skew Normal Mixtrue  ###
# Y: (n, p)-matrix of observations 
# K: number of components
# gam: value of gamma in the weight
# C: eigen constraint
# maxitr: maximum number of iteration 

RSNM <- function(Y, K=2, gam=0.2, C=10, maxitr=100, alpha=0.01, SE=T, init.fix=T){
  EPS <- 10^(-6)   # torelence rate
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  
  # skew normal density function 
  Density <- function(xx, Mu, Sig, Psi){
    Om <- array(NA,c(p, p, K))
    som <- array(NA,c(p, p, K))
    Alpha <- matrix(NA, p, K)
    for(k in 1:K){
      Om[,,k] <- Sig[,,k]+Psi[,k]%*%t(Psi[,k])
      som[,,k] <- diag( sqrt( diag(Om[,,k]) ) )
      denom <- as.vector( sqrt( 1-Psi[,k]%*%solve(Om[,,k])%*%Psi[,k] ) )
      Alpha[,k] <- as.vector( som[,,k]%*%solve(Om[,,k])%*%Psi[,k] )/denom
    }
    nn <- dim(xx)[1]
    f <- matrix(NA, nn, K)
    for(k in 1:K){
      f[,k] <- 2*dmvnorm(xx, Mu[,k], Om[,,k])*pnorm( Alpha[,k]%*%solve(som[,,k])%*%(t(xx)-Mu[,k]) )
    }
    return(f)
  }
  
  # Initial values (using "otrimle")
  if(init.fix){ set.seed(1) }
  PP <- rep(1/K, K)
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
  Para <- NA
  for(j in 1:maxitr){
    # save previous values
    val0 <- val
    Para0 <- Para
    Mu0 <- Mu
    
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
      ratio <- M[,k]/TT[,k]
      ratio[ratio<(-20)] <- -20
      a1 <- det(Sig[,,k])^(-gam/2)*(2*pi)^(-gam*p/2)*pnorm(ratio)/pnorm(ratio)*(TT[,k]/Tau[,k])
      a2 <- exp(-0.5*gam*c3 - 0.5*Delta[,k]^2/Tau[,k]^2 + 0.5*ratio^2 )
      U[,k] <- a1*a2
    }
    Dens <- Density(Y, Mu, Sig, Psi)   # density values
    uu <- t(t(Dens)*PP)
    uu <- uu/apply(uu,1,sum)   # posteiror probability
    
    # calculation of V
    V0 <- U
    V1 <- matrix(NA,N,K)
    V2 <- matrix(NA,N,K)
    for(k in 1:K){
      ratio <- M[,k]/TT[,k]
      #ratio[ratio<(-20)] <- -20
      Mill <- exp(dnorm(ratio, log=T) - pnorm(ratio, log.p=T))
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
      #num <- apply( uu[,k]*V1[,k]*t(t(Y)-Mu0[,k]), 2,sum)
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
    m.ast <- optim(par=1, fn=opt.M, method="L-BFGS-B", lower=0.001)$par
    
    # update Sigma
    for(k in 1:K){
      Sig[,,k] <- t(H[[k]])%*%diag(Trunc(D[,k],m.ast))%*%H[[k]]
      Sig[,,k] <- 0.5*(Sig[,,k]+t(Sig[,,k]))
    }
    
    ## making vector of parameters
    Para <- c()
    for(k in 1:K){
      est <- c(Mu[,k], Psi[,k], Sig[,,k][lower.tri(Sig[,,k])==F], PP[k])
      lab <- paste0("(", k, ")-", c(paste0("mu", 1:p), paste0("psi", 1:p), paste0("sig", 1:(p*(p+1)/2)), "mix"))
      names(est) <- lab 
      Para <- c(Para, est)
    }
    
    ## convergence check 
    Dens <- Density(Y, Mu, Sig, Psi)
    val <- mean( log( apply((t(Dens)*PP),2,sum) ))
    if(is.na(val)){ val <- 0 }
    if(abs(val)==Inf){ val <- 1000 }
    if( abs(val-val0)/(abs(val0)+0.01) < EPS ){ break }
  }

  ## function of estimating equations
  qq <- 2*p + p*(p+1)/2 + 1
  EE <- function(para){
    para <- matrix(para, qq, K)
    Mu <- para[1:p,]
    Psi <- para[(p+1):(2*p),]
    PP <- para[qq,]
    ss <- para[-c(1:(2*p), qq), ]
    Sig <- array(NA, c(p, p, K))
    for(k in 1:K){
      mat <- matrix(0, p, p)
      mat[lower.tri(mat)==F] <- ss[,k]
      mat <- mat + t(mat)
      diag(mat) <- diag(mat)/2
      Sig[,,k] <- mat
    }
  
    BK <- det(Sig[,,K])^(-gam/2)*(2*pi)^(-p*gam/2)*(1+gam)^(-p/2)
    E.mat <- c()
    for(k in 1:K){
      EE.mu <- matrix(NA, N, p)
      EE.sig <- matrix(NA, N, p*(p+1)/2)
      EE.psi <- matrix(NA, N, p)
      EE.pp <- c()
      gg <- gam*det(Sig[,,k])^(-gam/2)*(2*pi)^(-p*gam/2)*(1+gam)^(-p/2-1)
      B <- det(Sig[,,k])^(-gam/2)*(2*pi)^(-p*gam/2)*(1+gam)^(-p/2)
      for(i in 1:N){
        EE.mu[i, ] <- uu[i,k]*( V0[i,k]*(Y[i,]-Mu[,k]) - V1[i,k]*Psi[,k] )
        t1 <- uu[i,k]*V0[i,k]*Sig[,,k]
        t2 <- -uu[i,k]*V0[i,k]*(Y[i,]-Mu[,k])%*%t(Y[i,]-Mu[,k])
        t3 <- 2*uu[i,k]*V1[i,k]*(Y[i,]-Mu[,k])%*%t(Psi[,k])
        t4 <- -uu[i,k]*V2[i,k]*Psi[,k]%*%t(Psi[,k])
        t5 <- -gg*uu[i,k]*Sig[,,k]
        mat <- t1 + t2 + t3 + t4 + t5
        EE.sig[i,] <- mat[lower.tri(mat)==F]
        EE.psi[i,] <- uu[i,k]*( V1[i,k]*(Y[i,]-Mu[,k]) - V2[i,k]*Psi[,k] )
        EE.pp[i] <- uu[i,k]*V0[i,k]/PP[k]/B - uu[i,K]*V0[i,K]/PP[K]/BK
      }
      val <- cbind(EE.mu, EE.psi, EE.sig, EE.pp)
      dimnames(val)[[2]] <- lab 
      E.mat <- cbind(E.mat, val)
    }
    E.mat <- E.mat[,-qq*K]
    return(E.mat)
  }
  
  ##  varaince-covariance matrix
  if(SE==T){
    EE.c <- EE(Para)
    VG <- t(EE.c)%*%EE.c/N
    
    r <- length(Para)-1
    dG <- matrix(NA, r, r)
    for(j in 1:r){
      Para2 <- Para
      Para2[j] <- Para0[j]
      dG[j,] <- apply(EE(Para) - EE(Para2), 2, mean) / (Para[j] - Para0[j])
    }
    
    AV <- solve(dG)%*%VG%*%solve(t(dG))/N
    SE <- sqrt(diag(AV))
    names(SE) <- names(Para)[1:r]
  }
  
  ## outlier detection
  group <- apply(uu, 1, which.max)
  out <- list()
  for(k in 1:K){
    ID <- (1:N)[group==k]
    Nk <- length(ID)
    if(Nk>0){
      dens <- Density(Y, Mu, Sig, Psi)[group==k, k]
      MC <- 10000
      ref.rn <- matrix(NA, MC, p)
      for(i in 1:MC){
        ref.rn[i,] <- mvrnorm(1, Mu[,k], Sig[,,k]) + Psi[,k]*rtruncnorm(1, a=0)
      }
      ref.val <- Density(ref.rn, Mu, Sig, Psi)[,k]
      th <- quantile(ref.val, prob=alpha)
      out[[k]] <- ID[dens<th]
    }
  }
  out <- unlist(out)
  if(gam>0){ group[out] <- 0 }
  
  # BIC
  Dens <- Density(Y, Mu, Sig, Psi)
  dens <- apply(t(Dens)*PP, 2, sum)    # estimated density 
  ML <- N*mean( log(dens[-out]) )
  BIC <- -2*ML + log(N)*(2*p*K+K*p*(p+1)/2+K-1)
  
  # Summary
  Res <- list(Mu=Mu, Sigma=Sig, Psi=Psi, Pi=PP, out=out, BIC=BIC, cluster=group, post.prob=uu, itr=j, SE=SE)
  return(Res)
}


