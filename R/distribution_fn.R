#' log-likelihood function for problem 4.2
#' @param theta evaluated at theta = (alpha,beta,mu,lambda)
#' @param data observed counts
#' @export
em.loglike <- function(theta,counts)
{
  s <- 0
  for (i in 0:16)
  {
    s <- s + counts[i+1] * log(em.Pi(i,theta))
  }
  s
}

em.Pi <- function(i,theta)
{
  (i==0)*theta[1] + theta[2]*theta[3]^i*exp(-theta[3]) + (1-theta[1]-theta[2])*theta[4]^i*exp(-theta[4])
}

em.z0 <- function(theta)
{
  theta[1] / em.Pi(0,theta)
}

em.t <- function(i,theta)
{
  theta[2] * theta[3]^i * exp(-theta[3]) / em.Pi(i,theta)
}

em.p <- function(i,theta)
{
  (1-theta[1] - theta[2]) * theta[4]^i * exp(-theta[4]) / em.Pi(i,theta)
}
