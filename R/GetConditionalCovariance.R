#' get conditional covariance of one matrix given another:
#' @param mat Matrix for which conditional covariance is desired
#' @param condmat Matrix to be conditioned on
#' @param outputtype Either "cor" or "cov"
#' @return A conditional covariance (or correlation) matrix
#' @importFrom Rfast cova
get_conditional_correlation <- function(mat = NULL, condmat = NULL, outputtype = "cor") {
  cov_aa <- Rfast::cova(mat)
  cov_bb <- cov(condmat)
  cov_ab <- cov(mat, condmat)

  condcov <- cov_aa - cov_ab %*% solve(cov_bb) %*% t(cov_ab)

  if (outputtype == "cor") {
    condcov <- cov2cor(condcov)
  }
  return(condcov)
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
