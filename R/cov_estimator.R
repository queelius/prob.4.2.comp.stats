#' Bootstrap covariance estimator of EM point estimator
#' @param theta.em An EM point estimator of theta given observed counts
#' @param counts observed count data (n0,n1,...,n16)
#' @param m bootstrap replicates
#' @param debug whether to print out debugging info while running
#' @keywords EM algorithm
#' @export
em.cov_estimator.bootstrap <- function(theta.em, counts, m=10000, debug=T)
{
  #source("em_estimator.R")
  #source("data_convert.R")

  thetas <- rbind(theta.em)
  N <- sum(counts)
  data <- em.counts_to_responses(counts)
  for (i in 2:m)
  {
    indices <- sample(N,N,replace=T)
    thetas <- rbind(thetas,em(theta.em,em.responses_to_counts(data[indices]),eps,F))
    if (debug == T && i %% 500 == 0)
    {
      cat("iteration ", i, "\n")
      print(cov(thetas))
    }
  }
}

#' Observed information of EM point estimator based on an observed sample
#' @param theta.em An EM point estimator of theta given observed counts
#' @param data Observed sample of responses
#' @keywords observed information
#' @export
em.observed_info <- function(theta.em,data)
{
  #source("distribution_fn.R")
  library(numDeriv)
  l <- function(theta)
  {
    em.loglike(theta,data)
  }
  -hessian(l,theta.em)
}


#' Covariance matrix of EM point estimator based on the observed information matrix
#' @param theta.em An EM point estimator of theta given observed responses
#' @param data Observed sample of responses
#' @keywords covariance estimator
#' @export
em.cov_estimator.observed_info <- function(theta.em,data)
{
  #source("distribution_fn.R")
  library(numDeriv)
  solve(em.observed_info(theta.em,data))
}
