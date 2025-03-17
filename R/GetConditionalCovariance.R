#' get conditional covariance of one matrix given another:
#' @param mat Matrix for which conditional covariance is desired
#' @param condmat Matrix to be conditioned on
#' @param outputtype Either "cor" or "cov"
#' @return A conditional covariance (or correlation) matrix
get_conditional_correlation <- function(mat = NULL, condmat = NULL, outputtype = "cor") {
  cov_aa <- cova2(mat)
  cov_bb <- cov(condmat)
  cov_ab <- cov(mat, condmat)

  condcov <- cov_aa - cov_ab %*% solve(cov_bb) %*% t(cov_ab)

  if (outputtype == "cor") {
    condcov <- cov2cor(condcov)
  }
  return(condcov)
}

#' debugged Rfast::cova:
#' @param x matrix
#' @param center whether to center
#' @param large Same as Rfast::cova large argument
#' @importFrom Rfast colmeans
#' @importFrom Rfast Crossprod
#' @importFrom Rfast eachrow
#' @importFrom base crossprod
cova2 <- function(x, center = FALSE, large = FALSE) 
{
  n <- dim(x)[1]
  if (!center) {
    m <- sqrt(n) * Rfast::colmeans(x)
    if (large) {
      s <- (Rfast::Crossprod(x, x) - tcrossprod(m))/(n - 1)
    }
    else s <- (base::crossprod(x) - tcrossprod(m))/(n - 1)
  }
  else {
    m <- Rfast::colmeans(x)
    x <- Rfast::eachrow(x, m, oper = "-")
    if (large) {
      s <- Rfast::Crossprod(x, x)/(n - 1)
    }
    else s <- base::crossprod(x)/(n - 1)
  }
  s
}

#' convert covariance matrix to correlation:
#' @param covmat A covariance matrix
#' @return A correlation matrix
cov2cor <- function(covmat) {
  vars <- diag(covmat)
  cormat <- diag(vars^-0.5) %*% covmat %*% diag(vars^-0.5)
  dimnames(cormat) <- dimnames(covmat)
  return(cormat)
}
