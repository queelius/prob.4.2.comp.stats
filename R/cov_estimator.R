#' Bootstrap covariance estimator of EM point estimator
#'
#' Estimate the covariance of the EM point estimator for
#'     theta = (alpha,beta,mu,lamda)'
#' using Bootstrapping.
#' @param theta.em An EM point estimator of theta given observed counts
#' @param counts observed count data (n0,n1,...,n16)
#' @param m maximum bootstrap replicates
#' @param eps EM algorithm epsilon stopping condition
#' @param debug whether to print out debugging info while running
#' @export
em.cov.bs <- function(theta.em, counts, m=2000, eps=1e-6, debug=F)
{
  thetas <- rbind(theta.em)
  N <- sum(counts)
  data <- em.counts_to_responses(counts)
  for (i in 2:m)
  {
    indices <- sample(N,N,replace=T)
    resamp_counts <- em.responses_to_counts(data[indices])
    theta.bs <- em.estimator(theta.em,resamp_counts,eps,F)$estimate
    thetas <- rbind(thetas,theta.bs)
    if (debug == T && i %% 500 == 0)
    {
        cat("iteration", i, ": sample covariance:\n")
        print(cov(thetas))
    }
  }
  res <- cov(thetas)
  rownames(res) <- c("alpha","beta","mu","lamda")
  colnames(res) <- rownames(res)
  res
}

#' Standard error of EM point estimator based on Bootstrapping sample covariance
#' @param theta.em An EM point estimator of theta given observed responses
#' @param count Observed sample of counts
#' @export
em.sd.bs <- function(theta.em, counts, m=10000, em.eps=1e-6,bs.eps=1e-4,debug=F)
{
  res <- sqrt(diag(em.cov.bs(theta.em,counts,m,em.eps,bs.eps,debug)))
  colnames(res) <- c("alpha","beta","mu","lamda")
  res
}

#' Standard error of EM point estimator based on the observed information matrix
#' @param theta.em An EM point estimator of theta given observed responses
#' @param count Observed sample of counts
#' @export
em.sd.info <- function(theta.em, counts)
{
  res <- sqrt(diag(em.cov.info(theta.em,counts)))
  colnames(res) <- c("alpha","beta","mu","lamda")
  res
}

#' Covariance matrix of EM point estimator based on the observed information matrix
#' @param theta.em An EM point estimator of theta given observed responses
#' @param count Observed sample of counts
#' @export
em.cov.info <- function(theta.em,counts)
{
  res <- matrix(
    c(0,0,0,0,
    0,em.beta.var(theta.em,counts),0,0,
    0,0,em.mu.var(theta.em,counts),
    0,0,0,em.lamda.var(theta.em,counts)),nrow=4,byrow=T)
  rownames(res) <- c("alpha","beta","mu","lamda")
  colnames(res) <- rownames(res)
  res
}

em.beta.var <- function(theta.em,counts)
{
  s <- 0
  for (i in 0:16)
  {
    s <- s + counts[i+1]*((theta[3]^i*exp(-theta[3])-theta[4]^i*exp(-theta[4]))/my.Pi(i,theta))^2
  }
  1/s
}

em.mu.var <- function(theta.em,counts)
{
  s <- 0
  for (i in 0:16)
  {
    s <- s + counts[i+1]*((theta[2]*exp(-theta[3])*i*theta[3]^(i-1)-theta[2]*theta[3]^i*exp(-theta[3]))/my.Pi(i,theta))^2
  }
  1/s
}

em.lamda.var <- function(theta.em,counts)
{
  s <- 0
  gam <- 1-theta[1]-theta[2]
  for (i in 0:16)
  {
    s <- s + counts[i+1]*((gam*theta[4]^(i-1)*exp(-theta[4]) - gam*theta[4]^i*exp(-theta[4]))/my.Pi(i,theta))^2
  }
  1/s
}

