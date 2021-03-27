##############################################################################################
# This script assumes that the working directory has been set to the location of cod-contour.R
# Use getwd() to check the working directory, and setwd() to modify it.
##############################################################################################

require(ks)
rep <- 1000
loci <- 23
locs <- 1:loci
n <- 136
cutoff <- 15

rates_n <- 134
rates <- rep(NA, rates_n)
rates[1:11] <- format(0:10/10)
rates[12] <- 1.25
rates[13:18] <- format(seq(1.5,4,0.5))
rates[19:34] <- formatC(5:20,digits=1,format="f")
rates[35:38] <- formatC(seq(25,40,5),digits=1,format="f")
rates[39:rates_n] <- formatC(seq(50,1000,10),digits=1,format="f")
lambda_n <- 101

x_lab=expression(paste("Tail statistic ", bar(italic(Z))[italic(J)]*"("*italic(n)*","*italic(L)*","*italic(K)*")"),sep="")
y_lab=expression(paste("Singleton statistic ", bar(italic(Z))[italic(I)]*"("*italic(n)*","*italic(L)*","*italic(K)*")"),sep="")

exp_growth <- matrix(rep(NA, rates_n * 2 * rep), ncol = 2)
for (i in 1:rates_n) {
  filename <- paste("./sim-output/exp_",rates[i],".txt",sep="")
  tmp <- as.matrix(read.table(filename))[,1:cutoff]
  tmp[,2] <- 1 - rowSums(tmp[,1:(cutoff - 1)])
  for (j in 1:rep) {
    exp_growth[(i - 1) * rep + j,] <- colMeans(tmp[(j - 1) * loci + locs,2:1])
  }
}
lambda <- matrix(rep(NA, lambda_n * 2 * rep), ncol = 2)
for (i in 1:lambda_n) {
  filename <- paste("./sim-output/xi_", sprintf(1 + (i - 1) / (lambda_n - 1), fmt = '%#.2f'), ".txt", sep="")
  tmp <- as.matrix(read.table(filename))[,1:cutoff]
  tmp[,2] <- 1 - rowSums(tmp[,1:(cutoff - 1)])
  for (j in 1:rep) {
    lambda[(i - 1) * rep + j,] <- colMeans(tmp[(j - 1) * loci + locs,2:1])
  }
}

kde_exp <- kde(exp_growth)
kde_lambda <- kde(lambda)

pdf(paste("demo-contour.pdf",sep=""))
par(mar=c(5,5,4,2))
plot(kde_exp, cont = 99, xlim = c(0,1), ylim= c(0,1), xlab=x_lab,ylab=y_lab,drawlabels=FALSE,col="violet",las=1,cex.axis=1.5, cex.lab=2)
plot(kde_lambda, cont = c(25,50,75,99), add=TRUE,xaxt="n",yaxt="n",las=1,cex.axis=1.5, cex.lab=2,labcex=1,method="edge")

filename <- "./demo.sfs.out"
data <- as.vector(read.table(filename)[,2:n])
data <- as.numeric(data) / sum(data)
data <- c(data[1],sum(data[cutoff:(n-1)]))
points(data[2],data[1],col="green")
legend("topright",legend=c("Observed",expression("Beta("*2-alpha*","*alpha*")"),"Kingman"),pch=c(1,NA,NA),lty=c(NA,1,1),col=c("green","black","violet"),cex=1.7,lwd=c(1,3,3))
dev.off()