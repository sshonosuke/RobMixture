library(MASS)


logistic <- function(x){ 
  val <- 1/(1+exp(-x)) 
  val[val<10^(-8)] <- 10^(-8)
  val[val>0.99999] <- 0.99999
  return(val)
}

###  Robust Mixture-of-experts (linear regression)  ###
# K=2
# single covariate
RMOE <- function(Y, X, gam=0.2, maxitr=100, beta.init=NULL, num.init=10, alpha=0.01){
  K <- 2
  EPS <- 10^(-6)   # torelence rate
  N <- length(Y)
  XX <- cbind(1, X)

  # initial values of Beta
  Beta.init <- list()
  Beta.init[[1]] <- beta.init
  if( is.null(beta.init) ){
    for(l in 1:num.init){
      ind <- sort(sample(1:N, round(N/2)))
      Beta.init[[l]] <- cbind(coef(rlm(Y[ind]~X[ind])), coef(rlm(Y[-ind]~X[-ind])))
    }
  }
  
  LL <- length(Beta.init)  
  
 ### EEE algorithm
  Eta.set <- list()
  Beta.set <- list()
  Sig.set <- list()
  Out.set <- list()
  group.set <- list()
  BIC <- c()
  for(l in 1:LL){
    Beta <- Beta.init[[l]]
    Eta <- c(0, 0)
    Sig <- rep(sd(Y), 2)
    
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
      
      B <- c()
      G <- c()
      for(k in 1:K){
        B[k] <- (Sig[k]^2)^(-gam/2)*(2*pi)^(-gam/2)*(1+gam)^(-1/2)
        G[k] <- gam*B[k]/(1+gam)
      }
      
      ## EE-step
      opt <- function(eta){
        pp <- logistic(eta[1]+eta[2]*X)
        val <- -sum( t( t(uu*ww)/B )* cbind(log(pp), log(1-pp)) )
        val <- min(val, 10^10)
        return(val)
      }
      Eta <- optim(par=Eta, lower=c(-100,-100), upper=c(100,100), fn=opt, method="L-BFGS-B")$par
      
      for(k in 1:K){
        if(sum(uu[,k])>3){
          # Mean
          Beta[,k] <- solve(t(XX)%*%(uu[,k]*ww[,k]*XX))%*%t(XX)%*%(uu[,k]*ww[,k]*Y)
          # Variance 
          mu <- Beta[1,k]+Beta[2,k]*X
          SS <- sum(uu[,k]*ww[,k]*(Y-mu)^2) / ( sum(uu[,k]*ww[,k])-sum(uu[,k])*G[k] )
          if(SS<0){
            SS <- sum(uu[,k]*ww[,k]*(Y-mu)^2) / sum(uu[,k]*ww[,k])
          }
          Sig[k] <- sqrt(SS)
        }
      }
      
      # Convergence check
      dd <- sum(abs(Beta-Beta0)) / (sum(abs(Beta0))+0.001)
      if( dd < EPS ){ break }
    }
    Eta.set[[l]] <- Eta
    Sig.set[[l]] <- Sig
    Beta.set[[l]] <- Beta
    
    # outlier detection
    MC <- 10000
    group <- apply(uu, 1, which.max)
    out <- list()
    for(k in 1:K){
      ID <- (1:N)[group==k]
      Nk <- length(ID)
      if(Nk>0){
        mu <- as.vector(XX%*%Beta[,k])[group==k]
        dens <- dnorm(Y[group==k], mu, Sig[k])
        z <- Sig[k]*rnorm(MC)
        rn <- dnorm(z, 0, Sig[k])
        th <- quantile(rn, prob=alpha)
        out[[k]] <- ID[dens<th]
      }
    }
    out <- unlist(out)
    if(gam>0){ group[out] <- 0 }
    group.set[[l]] <- group
    Out.set[[l]] <- out
    
    # BIC
    Mu <- XX%*%Beta
    dens <- matrix(NA, N, K)
    for(k in 1:K){
      dens[,k] <- dnorm(Y, Mu[,k], Sig[k])
    }
    PP <- logistic(Eta[1] + Eta[2]*X)
    dens <- PP*dens[,1] + (1-PP)*dens[,2]    # estimated density 
    ML <- N*mean( log(dens[-out]) )
    BIC[l] <- -2*ML + log(N)*(2+K*3)
  }
  
  opt <- which.min(BIC)
  Beta <- Beta.set[[opt]]
  Eta <- Eta.set[[opt]]
  Sig <- Sig.set[[opt]]
  out <- Out.set[[opt]]
  group <- group.set[[opt]]
  
  # Summary
  Res <- list(beta=Beta, sigma=Sig, eta=Eta, out=out, group=group)
  return(Res)
}








