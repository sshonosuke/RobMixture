###-----------------------------------------------------------###
###        Simulation of Gaussian mixture                     ###
###-----------------------------------------------------------###
set.seed(123)
source("WCE-function.R")


N <- 500   # number of observations

K <- 3  # number of components
p <- 2  # dimension
om <- 0.05   # outlier ratio

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



# data generation
Y <- matrix(NA, N, p)
for(i in 1:N){
  ch <- sample(1:K, 1, prob=PP0)
  Y[i,] <- mvrnorm(1, Mu0[,ch], Sig0[,,ch])
}

# outlier 
for(i in (nb+1):N){
  dd <- 1
  while(max(dd)>0.005){
    prop <- c(runif(1,-10,10), runif(1,-5,5), runif(p-2,-3,3))
    Y[i,] <- prop
    dd <- rep()
    for(k in 1:3){  dd[k] <- dmvnorm(prop, Mu0[,k], Sig0[,,k])  }
  }
}

plot(Y, main="simulated data")





##  Estimation  ##
opt.K <- RGM.selection(Y, gam=0.5, C=20)
fit <- RGM(Y, K=opt.K, gam=0.5, C=20)

fit$Mu
fit$Sigma
fit$Pi


