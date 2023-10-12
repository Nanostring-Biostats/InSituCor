

#' Derive and score modules of spatially co-expressed genes
#'
#' Runs the complete workflow for deriving and scoring spatially co-expressed gene modules.
#' Namely: 1. computes the matrix of each cell's gene expression "environment",
#' 2. derives a matrix of genes' conditional correlation in the environment matrix, conditional
#' on factors like cell type and tissue ID,
#' 3. defines modules of co-expressed genes from this conditional correlation matrix,
#' and scores each cell for each module, both in its single cell profile and in its environment.
#' 4. For each module, scores each cell type for its role in driving the contribution to each gene.
#' @param counts Single cell expression matrix. Normalizing the data to give every cell the same total expression is preferred.
#' @param conditionon Data frame of variables to be conditioned on when computing gene correlation.
#'   Rows correspond to the rows of counts and xy. At a minimum, it is recommended to include cell type
#'   and tissue ID. Including cell signal strength (total counts) and background (negmean) is also recommended.
#' @param celltype Vector of cell type assignments
#' @param neighbors Neighbor relationships, stored as a sparse matrix
#' @param xy Matrix of xy coordinates.
#' @param k k for k-nearest neighbor network building
#' @param radius Radius for neighbor network building
#' @param tissue Used for neighbor network building. Neighbors will only be considered for cell with the same tissue value.
#' @param roundcortozero Correlation values below this threshold will be stored as 0 to allow for a sparse matrix
#' @param min_module_size Modules smaller than this are discarded
#' @param max_module_size Modules bigger than than this are subclustered
#' @param gene_weighting_rule How to define modules' gene weights, absed on gene expression levels.
#'   One of "inverse_sqrt", "inverse", or "identity".
#' @param resolution Argument to igraph::cluster_leiden. Lower values produce bigger clusters. 
#' @param corthresh Only correlations above this value will go into the adjacency graph fed into leiden clustering
#' @param min_module_cor Only keep modules with average cor above this value.
#' @param max_cells If there are more than this many cells, certain steps will use a random subset of this size.
#'   Output will still be for all cells.
#' @param attribution_subset_size Subsample size to use in attribution analysis. This function only extracts correlations, so 5000 cells is ample.
#' @param verbose Whether to print progress
#' @return A list, with the following elements:
#' \itemize{
#'  \item condcor: A sparse matrix holding genes' conditional correlations.
#'  \item modules: A data frame detailing gene module membership and weights
#'  \item scores_env: A matrix of cell * module environment scores
#'  \item scores_sc: A matrix of cell * module single cell scores
#'  \item attributionmats: A list of matrices holding attribution scores for each cell type * gene in each module.
#'  \item celltypeinvolvement: A matrix giving the maximum attribution score for each cell type in each module.
#' }
#' @export
#' @importFrom methods as
#' @importFrom stats cor cov lm model.matrix rnorm
#' @importFrom utils head
insitucor <- function(counts, conditionon = NULL, celltype,
                       neighbors = NULL, xy = NULL, k = NULL, radius = NULL, tissue = NULL, # args for neighbor definition
                       min_module_size = 3, max_module_size = 20,                           # args for module definition
                       resolution = 0.02, corthresh = 0.1, min_module_cor = 0.1, gene_weighting_rule = "inverse_sqrt",   # more args for module definition
                       roundcortozero = 0.1, max_cells = 5000,                              # args for controlling memory and compute
                       attribution_subset_size = 5000,                                      # args for cell type attribution scoring
                       verbose = TRUE) {

  ## get conditional covariance matrix
  temp <- calcSpatialCor (counts = counts, conditionon = conditionon,
                     neighbors = neighbors, xy = xy, k = k, radius = radius,
                     roundcortozero = roundcortozero, max_cells = max_cells)
  condcor <- temp$condcor
  env <- temp$env
  neighbors <- temp$neighbors
  rm(temp)

  ## derive modules
  if (verbose) {
    print("Defining modules")
  }
  modules <- defineModules(condcor = condcor,
                              env = env,
                              min_module_size = min_module_size,
                              max_module_size = max_module_size,
                              resolution = resolution, 
                              corthresh = corthresh,
                              min_module_cor = min_module_cor,
                              gene_weighting_rule = "inverse_sqrt")

  ## score modules
  if (verbose) {
    print("Calculating module scores for single cells and for cell neighborhoods")
  }
  temp <- scoreModules(counts = counts,
                             weights = modules$weights,
                             neighbors = neighbors)
  scores_env <- temp$scores_env
  scores_sc <- temp$scores_sc
  rm(temp)

  ## calculate cell-gene attribution scores
  if (verbose) {
    print("Running cell type attribution analysis")
  }
  attribres <- cellTypeAttribution(modulescores = scores_env,
                                                 weights = modules$weights,
                                                 counts = counts,
                                                 celltype = celltype,
                                                 neighbors = neighbors,
                                                 nsub = attribution_subset_size,
                                                 verbose = verbose)

  ## return everything useful
  out = list(modules = modules$weightsdf,
             scores_env = scores_env,
             scores_sc = scores_sc,
             condcor = condcor,
             attributionmats = attribres$attributionmats,
             celltypeinvolvement = attribres$involvescores)
  return(out)
}


