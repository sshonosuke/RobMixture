library(mclust)

## load function
source("RSNM-function.R")

## load data
data(banknote)

Y <- banknote[,-1]
Y <- apply(Y, 2, scale)
pairs(Y)
p <- dim(Y)[2]
qq <- 2*p + p*(p+1)/2 + 1


## monitoring BIC 
gam.set <- seq(0.1, 0.5, by=0.05)   # candidate for gamma
L1 <- length(gam.set)
K.set <- c(2, 3, 4)   # candidate for K
L2 <- length(K.set)

IC <- matrix(NA, L1, L2)
C <- 50
alpha <- 0.01    # probability threshold for outlier detection 
for(k in 1:L1){
  for(l in 1:L2){
    IC[k, l] <- RSNM(Y, gam=gam.set[k], K=K.set[l], maxitr=200, C=C, alpha=alpha, SE=F)$BIC
  }
}

IC <- cbind(NA, IC)
matplot(IC, type="b", col=c(1,1,2,3), lty=1, xaxt="n", xlab=expression(gamma), ylab="BIC")
axis(1, at=c(1, 3, 5, 7, 9, 11), label=gam.set[c(1, 3, 5, 7, 9, 11)])



## Proposed method (robust skew normal mixture)
fit1 <- RSNM(Y, gam=0.3, K=2, maxitr=200, C=C, alpha=0.001, init.fix=F)
fit1$out
Psi <- fit1$Psi
SE.psi <- matrix(c(fit1$SE, NA), qq, 2)[(p+1):(2*p),]     # standard error of the skewness parameter
est <- cbind(Psi, SE.psi)
est[,c(1,3,2,4)]

pairs(Y, col=fit1$cluster+1, pch=c(16,2,4)[fit1$cluster+1])    # clustering result


## standard skew normal mixture
fit0 <- RSNM(Y, gam=0, K=2, maxitr=200, C=C)


## Classificaiton
true <- c()
true[banknote$Class=="genuine"] <- 1
true[banknote$Class=="counterfeit"] <- 2

pred0 <- apply(fit0$post.prob, 1, which.max)
pred1 <- apply(fit1$post.prob, 1, which.max)

sum(ifelse(pred0!=true, 1, 0))    # number of mis-classification of the standard SNM
sum(ifelse(pred1!=true, 1, 0))    # number of mis-classification of the robust SNM

