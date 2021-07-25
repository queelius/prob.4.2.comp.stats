#' probability density function (pdf) for problem 4.2
#' @param x density evaluated at x given index theta
#' @param theta index of indexed family of density functions, theta = (alpha,beta,mu,lambda)
#' @keywords density
#' @export
em.pdf <- function(x,theta)
{
  theta[1]*as.numeric(x==0) * theta[2]*dpois(x,theta[3]) * (1-theta[1]-theta[2])*dpois(x,theta[4])
}

#' log-likelihood function for problem 4.2
#' @param theta log-likelihood evaluated at theta = (alpha,beta,mu,lambda)
#' @param data observed response data
#' @keywords log-likehood
#' @export
em.loglike <- function(theta,data)
{
  #s <- 0
  #for (x in data) { s <- s + log(pdf(x,theta)) }
  #s
  sum(log(em.pdf(data,theta)))
}

