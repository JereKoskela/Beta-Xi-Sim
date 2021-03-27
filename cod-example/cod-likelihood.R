#################################################################################################
# This script assumes that the working directory has been set to the location of cod-likelihood.R
# Use getwd() to check the working directory, and setwd() to modify it.
#################################################################################################

require(ks)
rep <- 1000
loci <- 23
locs <- 1:loci
n <- 136
cutoff <- 15
lambda_size <- 101
lambda <- matrix(rep(NA, lambda_size * 2 * rep), ncol = 2)
tmp <- matrix(rep(NA, lambda_size * loci * rep * 15), ncol = 15)

for (i in 1:lambda_size) {
  filename <- paste("./sim-output/xi_", sprintf(1 + (i - 1) / (lambda_n - 1), fmt = '%#.2f'), ".txt", sep="")
  tmp <- as.matrix(read.table(filename))[,1:15]
  tmp[,15] <- 1 - rowSums(tmp[,1:14])
  for (j in 1:rep) {
    lambda[(i - 1) * rep + j,] <- colMeans(tmp[(j - 1) * loci + locs, c(1,15)])
  }
  assign(paste("kdes_",i,sep=""),kde(lambda[((i - 1) * rep + 1):(i * rep),]))
}

lik <- rep(NA, lambda_size)
filename <- "./demo.sfs.out"
data <- as.vector(read.table(filename)[,2:n])
data <- as.numeric(data) / sum(data)
data <- c(data[1],sum(data[cutoff:(n-1)]))
for (i in 1:lambda_size) {
  lik[i] <- predict(eval(parse(text = paste("kdes_",i,sep=""))), x = data)
}

pdf("demo-likelihood.pdf")
par(mar=c(5,5,4,2))
plot(1+0:(lambda_size-1)/(lambda_size-1), lik, type="b",pch=1, ylim=c(0,ceiling(max(lik))), xlab=expression(alpha), ylab="Likelihood",las=1,cex.axis=1.5, cex.lab=2)
#legend("topright",legend=c(outgroup,lik_mod),pch=c(rep(19,4),1,2),col=c(cols,rep("black",2)), cex=2)
dev.off()