###-----------------------------------------------------------###
###        Simulation of Gaussian mixture                     ###
###-----------------------------------------------------------###
set.seed(123)
source("RGM-function.R")

N <- 500   # number of observations
K <- 3  # number of components
p <- 5  # dimension
om <- 0.06   # outlier ratio

nb <- N*(1-om)



# true parameters
cc <- 3
mu1 <- c(cc, cc, rep(0,p-2))
mu2 <- c(cc, -cc, rep(0,p-2))
mu3 <- c(-cc, -cc, rep(0,p-2))
Mu0 <- cbind(mu1, mu2, mu3)

sig1 <- diag(rep(1,p))
sig2 <- diag(rep(1,p))
sig3 <- diag(rep(1,p))
sig1[1:2,1:2] <- matrix(c(2,0.3,0.3,1), 2, 2)
sig2[1:2,1:2] <- matrix(c(1,-0.3,-0.3,1), 2, 2)
sig3[1:2,1:2] <- matrix(c(1,0.3,0.3,2), 2, 2)
Sig0 <- array(NA,c(p,p,K))
Sig0[,,1] <- sig1
Sig0[,,2] <- sig2
Sig0[,,3] <- sig3

PP0 <- c(0.3,0.3,0.4)



###  data generation
Y <- matrix(NA, N, p)     # data 
clust <- c()    # true classification (0 denotes outlier)

# non-outlier
for(i in 1:N){
  ch <- sample(1:K, 1, prob=PP0)
  Y[i,] <- mvrnorm(1, Mu0[,ch], Sig0[,,ch])
  clust[i] <- ch
}

# outlier
if(om>0){
  for(i in (nb+1):N){
    clust[i] <- 0
    dd <- 0
    while(dd < (5*p)){
      prop <- c(runif(1,-10,10), runif(1,-5,5), runif(p-2,-3,3))
      Y[i,] <- prop
      D <- rep()
      for(k in 1:3){
        D[k] <- t(prop-Mu0[,k])%*%solve(Sig0[,,k])%*%(prop-Mu0[,k])
      }
      dd <- min(D)
    }
  }
}


plot(Y, col=clust+1, main="simulated data")    # black points are outliers 




##  Estimation  ##
select <- RGM.select(Y, K.min=2, K.max=7, gam=0.2)
select$opt.K    # selected number of components

fit <- select$opt.fit
fit$Mu     # mean vector
fit$Sigma   # variance-covariance matrix
fit$Pi    # mixing proportion 
fit$out    # detected outliers


