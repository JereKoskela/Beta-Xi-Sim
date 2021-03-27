###################################################################################################
# This script assumes that the working directory has been set to the location of singleton-errors.R
# Use getwd() to check the working directory, and setwd() to modify it.
###################################################################################################

sample_sfs <- function(sfs, m, p, q, r) {
  n <- length(sfs) + 1
  s <- sum(sfs[2:(n-2)])
  A <- matrix(c(1-p-q-r, r, r-q, 1-r), nrow=2)
  A_inv <- solve(A)
  mu <- A_inv %*% (sfs[c(1,n-1)] - c((m-s)*q, 0))
  gamma <- rep(NA, 2)
  if (p * (1 - p) + r * (1 - r) < q * (1 - q)) {
    gamma[1] <- (m - s) * q * (1 - q)
  } else {
    gamma[1] <- (m - s) * (p * (1 - p) + r * (1 - r))
  }
  gamma[2] <- (m - s) * r * (1 - r)
  tmp <- gamma[1] * (1 - r)^2 + max(0, (r - q) * (2 - r - q) * gamma[2])
  gamma[2] <- r^2 * gamma[1] + (1 - p - q - r) * (1 - p - q + r) * gamma[2]
  gamma[1] <- tmp
  gamma <- gamma / (1 - p - q - r * (2 - p - 2 * q))^2
  M <- 1 / sqrt(r * (1 - r) * (m - s - 1) * q * (1 - q))
  sig <- rep(NA, 2)
  repeat {
    ret <- rnorm(2, mu, sqrt(gamma))
    sig[1] <- ret[1] * (p * (1 - p) + r * (1 - r) - q * (1 - q)) + ret[2] * (r * (1 - r) - q * (1 - q)) + (m - s) * q * (1 - q)
    sig[2] <- sum(ret) * r * (1 - r)
    covar <- matrix(c(sig[1], -sig[2], -sig[2], sig[2]), nrow=2)
    alpha <- exp(-t(ret - mu) %*% (t(A) %*% solve(covar) %*% A - solve(diag(gamma))) %*% (ret - mu) / 2 - log(det(covar)) / 2 - log(M))
    if (runif(1) < alpha) {
      break
    }
  }
  ret
}

rep <- 100
filename <- "../demo.sfs.out"
cod_sfs <- as.vector(read.table(filename))
n <- length(cod_sfs) - 1
simulated_error_sfs <- matrix(rep(NA, rep * (n - 1)), nrow=rep)
m <- sum(cod_sfs)
cod_sfs <- as.numeric(cod_sfs[2:n])
for (i in 1:rep) {
  simulated_error_sfs[i,] <- cod_sfs
}
###########################################################################################
# Desired error probabilities. The script will loop over all combinations and produce one
# plot for each.

# Probabilities of missing a true singleton
ps <- c(1e-3, 1e-2, 1e-1)

# Probabilities of calling a false singleton
qs <- c(1e-4, 1e-3, 1e-2)

# Probabilities of an outgroup error resulting in a singleton called as an anti-singleton,
# or vice versa
rs <- c(1e-6, 1e-5)
###########################################################################################
for (p in ps) {
  for (q in qs) {
    for (r in rs) {
      for (i in 1:rep) {
        simulated_error_sfs[i,] <- cod_sfs
        simulated_error_sfs[i,c(1,n-1)] <- sample_sfs(cod_sfs, m, p, q, r)
      }
      simulated_error_sfs <- simulated_error_sfs / rowSums(simulated_error_sfs)

      pdf(paste("sfs-err-p=",p,"-q=",q,"-r=",r,".pdf", sep=""))
      par(mar=c(5,5,4,2))
      plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), log(simulated_error_sfs[1,]) - log(1 - simulated_error_sfs[1,]), type="l",ylim=c(-9,1),las=1,cex.axis=1.5, main=paste0("p = ", p, ", q = ", q, ", r = ", r), cex.main=2, cex.lab=2, xlab="Logit allele frequency", ylab="Logit normalised SFS")
      par(new=TRUE)
      for (i in 2:rep) {
        plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), log(simulated_error_sfs[i,]) - log(1 - simulated_error_sfs[i,]), type="l",ylim=c(-9,1),xlab="",ylab="",xaxt="n",yaxt="n")
        par(new=TRUE)
      }
      plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), log(cod_sfs/ sum(cod_sfs)) - log(1 - cod_sfs/ sum(cod_sfs)), type="l", col="red",ylim=c(-9,1),xlab="",ylab="",xaxt="n",yaxt="n")
      legend("topright", legend=c("observed", "simulated errors"), lty=c(1,1), lwd=c(3,3), col=c("red", "black"), cex=1.7)
      fig_lab="z"
      if (p == 1e-3 && q == 1e-3 && r == 1e-5) {
        fig_lab="a"
      }
      if (p == 1e-3 && q == 1e-2 && r == 1e-6) {
        fig_lab="b"
      }
      if (p == 1e-3 && q == 1e-4 && r == 1e-5) {
        fig_lab="c"
      }
      if (p == 1e-2 && q == 1e-2 && r == 1e-5) {
        fig_lab="d"
      }
      if (p == 1e-1 && q == 1e-2 && r == 1e-5) {
        fig_lab="e"
      }
      if (p == 1e-1 && q == 1e-4 && r == 1e-6) {
        fig_lab="f"
      }
      mtext(text=fig_lab, side=3, adj=0, outer=FALSE, cex=3)
      dev.off()
    }
  }
}