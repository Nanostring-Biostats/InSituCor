% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuildConditionalMatrix.R
\name{build_conditional_matrix}
\alias{build_conditional_matrix}
\title{Combine all variables to be conditioned on into a viable model matrix}
\usage{
build_conditional_matrix(variables)
}
\arguments{
\item{variables}{A list containing all the variables (could be vectors or matrices) that you wish to condition on}
}
\value{
A numeric matrix of continuous and indicator variables, ready to be conditioned on.
}
\description{
Combine all variables to be conditioned on into a viable model matrix
}
