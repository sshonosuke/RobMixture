###-----------------------------------------------------------###
###        Simulation of skew normal mixture                  ###
###-----------------------------------------------------------###
set.seed(123)
source("WCE-function.R")
library(truncnorm)


# density of skew-normal 
dens <- function(x1, x2, Mu, Sig, Psi){
  xx <- c(x1,x2)
  Om <- Sig+Psi%*%t(Psi)
  som <- diag( sqrt( diag(Om) ) )
  denom <- as.vector( sqrt( 1-Psi%*%solve(Om)%*%Psi ) )
  Alpha <- as.vector( som%*%solve(Om)%*%Psi )/denom
  2*dmvnorm(xx,Mu,Om)*as.vector( pnorm( Alpha%*%solve(som)%*%(xx-Mu) ) )
}


# settings
N <- 500   # number of observations
om <- 0.06    # outlier ratio
ran <- 15      # range of outliers

nb <- N*(1-om)


# true parameters
mu <- c(1.5, 1.5)
Mu0 <- cbind(mu, -mu)
Sig0 <- matrix(c(1,0.3,0.3,1), 2, 2)
psi <- c(2, 2)
Psi0 <- cbind(psi,-psi)
PP0 <- c(0.4, 0.6)


# data generation 
Y <- matrix(NA, N, 2)
for(i in 1:N){
  sn <- rtruncnorm(1,a=0)
  ep <- mvrnorm(1,c(0,0),Sig0)
  if(PP0[1]>runif(1)){
    Y[i,] <- Mu0[,1]+Psi0[,1]*sn+ep
  } else{
    Y[i,] <- Mu0[,2]+Psi0[,2]*sn+ep
  }
}

# outlier generation 
for(i in (nb+1):N){
  d1 <- d2 <- 1
  while(d1>0.005 | d2>0.005){
    prop <- c(runif(1, -ran, ran), runif(1, -ran, ran))
    d1 <- dens(prop[1], prop[2], Mu0[,1], Sig0, Psi0[,1])
    d2 <- dens(prop[1], prop[2], Mu0[,2], Sig0, Psi0[,2])
  }
  Y[i,]=prop
}

plot(Y, main="simulated data")





##  Estimation  ##
fit1 <- RSNM(Y, K=2, gam=0.5, C=20)
Mu1 <- fit1$Mu
Sig1 <- fit1$Sigma
Psi1 <- fit1$Psi
PP1 <- fit1$Pi

fit2 <- RSNM(Y, K=2, gam=0, C=20)
Mu2 <- fit2$Mu
Sig2 <- fit2$Sigma
Psi2 <- fit2$Psi
PP2 <- fit2$Pi


# contour plot of estimated density
mc <- 200     # number of evaluation points
x1 <- x2 <- seq(-10, 10, length=mc)     # evaluation points

# estimated density (robust skew normal mixture)
Z1 <- matrix(NA, mc, mc)
Z2 <- matrix(NA, mc, mc)
for(i in 1:mc){
  for(j in 1:mc){ 
    Z1[i,j] <- dens(x1[i], x2[j], Mu1[,1], Sig1[,,1], Psi1[,1]) 
    Z2[i,j] <- dens(x1[i], x2[j], Mu1[,2], Sig1[,,2], Psi1[,2]) 
  }
}
Z.est1 <- PP1[1]*Z1 + PP1[2]*Z2


# estimated density (standard skew normal mixture)
Z1 <- matrix(NA, mc, mc)
Z2 <- matrix(NA, mc, mc)
for(i in 1:mc){
  for(j in 1:mc){ 
    Z1[i,j] <- dens(x1[i], x2[j], Mu2[,1], Sig2[,,1], Psi2[,1]) 
    Z2[i,j] <- dens(x1[i], x2[j], Mu2[,2], Sig2[,,2], Psi2[,2]) 
  }
}
Z.est2 <- PP2[1]*Z1 + PP2[2]*Z2

# true density
tZ1 <- matrix(NA, mc, mc)
tZ2 <- matrix(NA, mc, mc)
for(i in 1:mc){
  for(j in 1:mc){ 
    tZ1[i,j] <- dens(x1[i], x2[j], Mu0[,1], Sig0, Psi0[,1]) 
    tZ2[i,j] <- dens(x1[i], x2[j], Mu0[,2], Sig0, Psi0[,2]) 
  }
}
tZ <- PP0[1]*tZ1 + PP0[2]*tZ2


contour(x1, x2, Z.est1, col="blue", lwd=1.5, main="skew normal mixture")
contour(x1, x2, Z.est2, col="red", lwd=1.5, main="", add=T)
contour(x1, x2, tZ, col="black", lwd=1.5, add=T)
points(Y[1:nb,], col=rgb(0, 0, 0, alpha=0.3))
points(Y[(nb+1):N,], col=rgb(0, 0, 0, alpha=0.3), pch=4)
legend("topleft", legend=c("WCE", "SNM", "True"), lty=1, col=c("blue", "red", "Black"))
legend("bottomright", legend=c("non-outlier", "outlier"), col=1, pch=c(1,4))
