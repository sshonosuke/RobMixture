library(mixtools)

# load data
data(tonedata)
X <- tonedata[,1]
Y <- tonedata[,2]

# load function
source("WCE-function.R")


# fit the robust mixture of experts
fit1 <- RMOE(Y, X, gam=0.5)
fit1$Beta
fit1$Eta

# fit the standard mixture of experts
fit2 <- hmeEM(Y, X)
fit2$beta


# sensitivity analysis (by adding outliers) 
EX <- c(X, rep(0, 10))
EY <- c(Y, rep(4, 10))

efit1 <- RMOE(EY, EX, gam=0.5)
efit1$Beta

efit2 <- hmeEM(EY, EX)
efit2$beta



# plot 
par(mfcol=c(1,2))

plot(EX, EY, main="RMOE", xlab="x", ylab="y", ylim=c(1.3, 4.1))
abline(efit1$Beta[,1], col=1)
abline(efit1$Beta[,2], col=2)
legend("topright", legend=c("cluster 1", "cluster 2"), lty=1, col=c(1,2))

plot(EX, EY, main="MOE", xlab="x", ylab="y", ylim=c(1.3, 4.1))
abline(efit2$beta[,1], col=1)
abline(efit2$beta[,2], col=2)
legend("topright", legend=c("cluster 1", "cluster 2"), lty=1, col=c(1,2))



