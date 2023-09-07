

#' Combine all variables to be conditioned on into a viable model matrix
#' @param variables A list containing all the variables (could be vectors or matrices) that you wish to condition on
#' @return A numeric matrix of continuous and indicator variables, ready to be conditioned on.
build_conditional_matrix <- function(variables) {
  # get number of obs:
  if (is.matrix(variables[[1]])) {
    nobs <- nrow(variables[[1]])
  }
  if (is.vector(variables[[1]])) {
    nobs <- length(variables[[1]])
  }
  # use the built-in lm machinery to define a model matrix for each variable:
  tempy = rnorm(nobs)
  allmat <- NULL
  for (i in 1:length(variables)) {
    if (is.factor(variables[[i]])) {
      variables[[i]] <- as.character(variables[[i]])
    }
    mod <- lm(tempy ~ variables[[i]])
    newmat <- model.matrix(mod)[, -1, drop = FALSE]
    # if it's not positive definite, remove a row:
    if (min(eigen(cov(newmat))$values) < 1e-4) {
      newmat <- newmat[, -1]
    }
    allmat <- cbind(allmat, newmat)
  }
  colnames(allmat) <- gsub("variables\\[\\[i\\]\\]", "", colnames(allmat))
  return(allmat)
}
