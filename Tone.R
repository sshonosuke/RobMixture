library(mixtools)

## load data
data(tonedata)
X <- tonedata[,1]
Y <- tonedata[,2]

## load function
source("RMOE-function.R")


## fit the robust mixture of experts
fit1 <- RMOE(Y, X, gam=0.3)
fit1$beta
fit1$eta
fit1$out    # detected outliers
hg1 <- fit1$group


## fit the standard mixture of experts
fit2 <- hmeEM(Y, X)
fit2$beta
hg2 <- apply(fit2$posterior, 1, which.max)


## sensitivity analysis (by adding outliers) 
EX <- c(X, rep(0, 10))
EY <- c(Y, rep(4, 10))


# robust method
efit1 <- RMOE(EY, EX, gam=0.3, beta.init=fit1$beta, alpha=0.005)
efit1$beta
efit1$eta
efit1$out
ehg1 <- efit1$group

# standard method
efit2 <- hmeEM(EY, EX)
efit2$beta
ehg2 <- apply(efit2$posterior, 1, which.max)


## Figure 
par(mfcol=c(2,2))
plot(X, Y, main="RMOE", xlab="x", ylab="y", col=c(4,2,1)[hg1+1], pch=ifelse(hg1>0, 1, 4))
abline(fit1$beta[,1], col=2)
abline(fit1$beta[,2], col=1)
legend("topleft", legend=c("cluster 1", "cluster 2", "outlier"), lty=c(1,1,NA), col=c(2,1,4),pch=c(1,1,4))

plot(EX, EY, main="RMOE", xlab="x", ylab="y", ylim=c(1.3, 4.1), col=c(4,2,1)[ehg1+1], pch=ifelse(ehg1>0, 1, 4))
abline(efit1$beta[,1], col=2)
abline(efit1$beta[,2], col=1)
legend("topright", legend=c("cluster 1", "cluster 2", "outlier"), lty=c(1,1,NA), col=c(2,1,4),pch=c(1,1,4))

plot(X, Y, main="MOE", xlab="x", ylab="y", col=hg2)
abline(fit2$beta[,1], col=1)
abline(fit2$beta[,2], col=2)
legend("topleft", legend=c("cluster 1", "cluster 2"), lty=1, col=c(1,2), pch=1)

plot(EX, EY, main="MOE", xlab="x", ylab="y", ylim=c(1.3, 4.1), col=ehg2)
abline(efit2$beta[,1], col=1)
abline(efit2$beta[,2], col=2)
legend("topright", legend=c("cluster 1", "cluster 2"), lty=1, col=c(1,2), pch=1)

