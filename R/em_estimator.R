#' EM algorithm
#'
#' EM algorithm estimator for problem 4.2
#' @param theta a starting guess for theta = (alpha,beta,mu,lambda)
#' @param counts observed count data (n0,n1,...,n16)
#' @param eps stopping condition
#' @param debug whether to print out debugging info while running
#' @export
em.estimator <- function(theta,counts,eps=1e-6,debug=T)
{
  update <- function(theta,counts)
  {
    N <- sum(counts)
    alpha <- counts[1] * em.z0(theta) / N
    beta <- 0
    mu_num <- 0
    mu_denom <- 0

    lam_num <- 0
    lam_denom <- 0

    for (i in 0:16)
    {
      ti <- em.t(i,theta)
      pi <- em.p(i,theta)

      beta <- beta + counts[i+1] * ti

      mu_num <- mu_num + i * counts[i+1] * ti
      mu_denom <- mu_denom + counts[i+1] * ti

      lam_num <- lam_num + i * counts[i+1] * pi
      lam_denom <- lam_denom + counts[i+1] * pi
    }

    beta <- beta / N
    mu <- mu_num / mu_denom
    lam <- lam_num / lam_denom

    c(alpha,beta,mu,lam)
  }

  theta.new <- NULL
  i <- 0
  repeat
  {
    i <- i + 1
    theta.new <- update(theta,counts)
    if (debug==T && i %% 100 == 0) { cat("iteration =",i," theta = (",theta.new,")'\n") }
    if (max(abs(theta.new - theta)) < eps) { break }
    theta <- theta.new
  }
  list(estimate=theta.new,iterations=i)
}
