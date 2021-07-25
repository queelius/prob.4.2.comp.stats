#' Bootstrap covariance estimator of EM point estimator
#' @param theta.em An EM point estimator of theta given observed counts
#' @param counts observed count data (n0,n1,...,n16)
#' @param m maximum bootstrap replicates
#' @param em.eps EM algorithm epsilon stopping condition
#' @param bs.eps Bootstrap covariance epsilon stopping condition
#' @param debug whether to print out debugging info while running
#' @keywords EM algorithm
#' @export
em.cov.bs <- function(theta.em, counts, m=10000, em.eps=1e-6,bs.eps=1e-4,debug=F)
{
  thetas <- rbind(theta.em)
  N <- sum(counts)
  data <- em.counts_to_responses(counts)
  cov.old <- NULL
  k <- 500

  for (i in 2:m)
  {
    indices <- sample(N,N,replace=T)
    resamp_counts <- em.responses_to_counts(data[indices])
    theta.bs <- em.estimator(theta.em,resamp_counts,em.eps,F)$estimate
    thetas <- rbind(thetas,theta.bs)
    if (i %% k == 0)
    {
      cov.new <- cov(thetas)
      if (i > k)
      {
        if (debug == T)
        {
          cat("iteration", i, ": cov delta is:\n")
          print(cov.new - cov.old)
        }
        if (max(abs(cov.new - cov.old)) < bs.eps) { break }
      }
      cov.old <- cov.new
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
#' @keywords observed information
#' @export
em.observed_info <- function(theta.em,data)
{
  library(numDeriv)
  l <- function(theta) { em.loglike(theta,data) }
  -hessian(l,theta.em)
}

#' Covariance matrix of EM point estimator based on the observed information matrix
#' @param theta.em An EM point estimator of theta given observed responses
#' @param data Observed sample of responses
#' @keywords covariance estimator
#' @export
em.cov.info <- function(theta.em,counts)
{
  data <- em.counts_to_responses(counts)
  obs_info <- em.observed_info(theta.em,data)
  solve(obs_info)
}

