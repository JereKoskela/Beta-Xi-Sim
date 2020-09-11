###################################################################################################
# This script assumes that the working directory has been set to the location of singleton-errors.R
# Use getwd() to check the working directory, and setwd() to modify it.
###################################################################################################

sample_xi <- function(sfs, m, p, q) {
  n <- length(sfs) + 1
  mu <- (sfs[1] - (m - sum(sfs[2:(n-1)])) * q) / (1 - p - q)
  sd_proposal <- sqrt((m - sum(sfs[2:(n-1)])) * max(p * (1 - p), q * (1 - q))) / (1 - p - q)
  M <- sqrt(max(p * (1 - p) / (q * (1 - q)), q * (1 - q) / (p * (1 - p))))
  repeat {
    ret <- rnorm(1, mu, sd_proposal)
    sd_target <- sqrt(ret * (p * (1 - p) - q * (1 - q)) + (m - sum(sfs[2:(n-1)])) * q * (1 - q)) / (1 - p - q)
    if (p == 0 || q == 0) {
      alpha <- dnorm(ret, mu, sd_target) / (dnorm(ret, mu, sd_proposal) * sqrt(m - sum(sfs[2:(n-1)])))
    } else {
      alpha <- dnorm(ret, mu, sd_target) / (dnorm(ret, mu, sd_proposal) * M)
    }
    if (runif(1) < alpha) {
      break
    }
  }
  ret
}

###########################################################################################
# Change these two parameters and rerun the script to produce plots.

p <- 0.9 # per-site probability of missing a true singleton
q <- 0.01 # per-site probability of recording a false singleton from a non-segregating site

###########################################################################################

simulated_bs_sfs <- colMeans(read.table("../sim-output/xi_1.txt"))
n <- length(simulated_bs_sfs) + 1
rep <- 100
simulated_error_sfs <- matrix(rep(NA, rep * (n - 1)), nrow=rep)
filename <- "../demo.sfs.out"
cod_sfs <- as.vector(read.table(filename))
m <- sum(cod_sfs)
cod_sfs <- as.numeric(cod_sfs[2:n])
for (i in 1:rep) {
  simulated_error_sfs[i,] <- cod_sfs
}

for (i in 1:rep) {
  simulated_error_sfs[i,1] <- sample_xi(cod_sfs, m, p, q)
}
simulated_error_sfs <- simulated_error_sfs / rowSums(simulated_error_sfs)
cod_sfs <- cod_sfs / sum(cod_sfs)

pdf(paste("bolthausen-sznitman-err-p=",p,"-q=",q,".pdf", sep=""))
plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), log(simulated_bs_sfs) - log(1 - simulated_bs_sfs), type="l", col="red", ylim=c(-9,1),xlab="",ylab="",las=1,cex.axis=1.3)
par(new=TRUE)
for (i in 1:rep) {
  plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), log(simulated_error_sfs[i,]) - log(1 - simulated_error_sfs[i,]), type="l", ylim=c(-9,1),xlab="",ylab="",xaxt="n", yaxt="n")
  par(new=TRUE)
}
plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), log(cod_sfs) - log(1 - cod_sfs), type="l", ylim=c(-9,1),col="blue",xlab="",ylab="",xaxt="n", yaxt="n")
legend("topright", legend=c("simulated","observed", "singleton errors"), pch=rep(16,3), col=c("red","blue", "black"), cex=1.3)
dev.off()