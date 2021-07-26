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
  cov(thetas)
}

em.sd.bs <- function(theta.em, counts, m=10000, em.eps=1e-6,bs.eps=1e-4,debug=F)
{
  sqrt(diag(em.cov.bs(theta.em,counts,m,em.eps,bs.eps,debug)))
}

#' Observed information of EM point estimator based on an observed sample
#' @param theta.em An EM point estimator of theta given observed counts
#' @param data Observed sample of responses
#' @export
em.observed_info <- function(theta.em,counts)
{
  library(numDeriv)
  -hessian(function(theta) { em.loglike(theta,em.counts_to_responses(counts)) },theta.em)
}

#' Covariance matrix of EM point estimator based on the observed information matrix
#' @param theta.em An EM point estimator of theta given observed responses
#' @param data Observed sample of responses
#' @export
em.cov.info <- function(theta.em,counts)
{
  solve(em.observed_info(theta.em,counts))
}

