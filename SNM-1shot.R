###-----------------------------------------------------------###
###        Simulation of skew normal mixture                  ###
###-----------------------------------------------------------###
set.seed(123)
source("RSNM-function.R")
library(truncnorm)


## density of skew-normal 
dens <- function(x1, x2, Mu, Sig, Psi){
  xx <- c(x1,x2)
  Om <- Sig+Psi%*%t(Psi)
  som <- diag( sqrt( diag(Om) ) )
  denom <- as.vector( sqrt( 1-Psi%*%solve(Om)%*%Psi ) )
  Alpha <- as.vector( som%*%solve(Om)%*%Psi )/denom
  2*dmvnorm(xx,Mu,Om)*as.vector( pnorm( Alpha%*%solve(som)%*%(xx-Mu) ) )
}


## settings
N <- 500   # number of observations
om <- 0.06    # outlier ratio
ran <- 10   # outlier range
nb <- N*(1-om)
p <- 2    # dimension
K <- 2
cc <- 1


## true parameters
mu <- c(cc, cc, rep(0, p-2))
Mu0 <- cbind(mu, -mu)
Sig0 <- diag( rep(1,p) )
Sig0[1:2,1:2] <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
psi <- c(2, 2, rep(0, p-2))
Psi0 <- cbind(psi, -psi)
PP0 <- c(0.4, 0.6)


## data generation 
Y <- matrix(NA, N, p)
clust <- c()

# genuine observations
for(i in 1:N){
  sn <- rtruncnorm(1, a=0)
  ep <- mvrnorm(1, rep(0, p), Sig0)
  ch <- sample(1:K, size=1, prob=PP0)
  clust[i] <- ch
  Y[i,] <- Mu0[,ch]+Psi0[,ch]*sn+ep
}

# outlier generation 
for(i in (nb+1):N){
  d1 <- d2 <- 1
  while(d1<5*p | d2<5*p){
    prop <- c(runif(1, -ran, ran), runif(1, -ran, ran), runif(p-2, -3, 3))
    d1 <- t(prop-Mu0[,1])%*%solve(Sig0)%*%(prop-Mu0[,1])
    d2 <- t(prop-Mu0[,2])%*%solve(Sig0)%*%(prop-Mu0[,2])
  }
  Y[i,] <- prop
  clust[i] <- 0
}


par(mfcol=c(1,1))
plot(Y, main="simulated data")





##  Estimation  ##
fit1 <- RSNM(Y, K=2, gam=0.2, C=10)
Mu1 <- fit1$Mu
Sig1 <- fit1$Sigma
Psi1 <- fit1$Psi
PP1 <- fit1$Pi

fit2 <- RSNM(Y, K=2, gam=0, C=10)
Mu2 <- fit2$Mu
Sig2 <- fit2$Sigma
Psi2 <- fit2$Psi
PP2 <- fit2$Pi



###  contour plot of estimated density
mc <- 200     # number of evaluation points
x1 <- x2 <- seq(-8, 8, length=mc)     # evaluation points

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



## Figure
contour(x1, x2, Z.est1, col="red", lwd=1.5, main="skew normal mixture")
contour(x1, x2, Z.est2, col="green", lwd=1.5, main="", add=T)
contour(x1, x2, tZ, col="blue", lwd=1.5, add=T)
points(Y[1:nb,], col=rgb(0, 0, 0, alpha=0.3))
points(Y[(nb+1):N,], col=rgb(0, 0, 0, alpha=0.3), pch=4)
legend("topleft", legend=c("WCE", "SNM", "True"), lty=1, col=c("red", "green", "blue"))
legend("bottomright", legend=c("non-outlier", "outlier"), col=1, pch=c(1,4))
