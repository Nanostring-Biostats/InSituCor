#' Small example CosMx data from a kidney sample
#'
#' A 2-FOV, 11-gene excerpt from a 1000-plex CosMx run on a kidney biopsy.
#'
#' @format A list with the following elements:
#'  \itemize{
#'  \item counts A matrix of raw counts, with cells i rows and genes in columns
#'  \item annot A data frame holding cell metadata
#'  \item xy A matrix of cell xy positions
#'  }
"cosmx_kidney"

#' Simulated example CosMx data from a kidney sample
#'
#' Data is simulated so most genes have no spatial pattern (beyond what cell type explains), while a few are highly spatially correlated. 
#'
#' @format A list with the following elements:
#'  \itemize{
#'  \item counts A matrix of raw counts, with cells i rows and genes in columns
#'  \item annot A data frame holding cell metadata
#'  \item xy A matrix of cell xy positions
#'  \item celltype Vector of cell types
#'  }
"simdat"
