############################################################################################
# This script assumes that the working directory has been set to the location of cod-power.R
# Use getwd() to check the working directory, and setwd() to modify it.
############################################################################################

require(ks)
rep <- 1000
loci <- 23
locs <- 1:loci
n <- 136
rates_n <- 134
rates <- rep(NA, rates_n)
rates[1:11] <- format(0:10/10)
rates[12] <- 1.25
rates[13:18] <- format(seq(1.5,4,0.5))
rates[19:34] <- formatC(5:20,digits=1,format="f")
rates[35:38] <- formatC(seq(25,40,5),digits=1,format="f")
rates[39:rates_n] <- formatC(seq(50,1000,10),digits=1,format="f")
lambda_n <- 101

exp_growth <- matrix(rep(NA, rates_n * 2 * rep), ncol = 2)
for (i in 1:rates_n) {
  filename <- paste("./sim-output/exp_",rates[i],".txt",sep="")
  tmp <- as.matrix(read.table(filename))[,1:15]
  tmp[,2] <- 1 - rowSums(tmp[,1:14])
  for (j in 1:rep) {
    exp_growth[(i - 1) * rep + j,] <- colMeans(tmp[(j - 1) * loci + locs,1:2])
  }
}

lambda <- matrix(rep(NA, lambda_n * 2 * rep), ncol = 2)
for (i in 1:lambda_n) {
  filename <- paste("./sim-output/xi_", sprintf(1 + (i - 1) / (lambda_n - 1), fmt = '%#.2f'), ".txt", sep="")
  tmp <- as.matrix(read.table(filename))[,1:15]
  tmp[,2] <- 1 - rowSums(tmp[,1:14])
  for (j in 1:rep) {
    lambda[(i - 1) * rep + j,] <- colMeans(tmp[(j - 1) * loci + locs,1:2])
  }
}

exp_kde <- kde(exp_growth)
lambda_kde <- kde(lambda)
q <- 0
for (i in 1:rates_n) {
  q <- max(q, quantile(predict(lambda_kde, x=exp_growth[((i - 1) * rep + 1):(i * rep),]) / pmax(predict(exp_kde, x=exp_growth[((i - 1) * rep + 1):(i * rep),]), .Machine$double.xmin), 0.99))
}
power <- rep(NA, lambda_n)
for (i in 1:lambda_n) {
  power[i] <- sum(predict(lambda_kde, x=lambda[((i - 1) * rep + 1):(i * rep),]) / pmax(predict(exp_kde, x=lambda[((i - 1) * rep + 1):(i * rep),]), .Machine$double.xmin) > q) / rep
}

pdf("demo-power.pdf")
par(mar=c(5,5,4,2))
plot(1+0:(lambda_n-1)/(lambda_n-1), power,type="l",xlab=expression(alpha), ylab="Power",xlim=c(1,2),ylim=c(0,1), las=1, cex.axis=1.5,cex.lab=2)
dev.off()