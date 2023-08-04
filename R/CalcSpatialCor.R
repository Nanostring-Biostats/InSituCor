
# sub-wrapper 1: from raw inputs to a conditional cor matrix:

#' Get conditional correlation matrix
#'
#' Builds a neighbor network, calculates an environment matrix, calculates and environment
#'  matrix of variables to be conditioned on, then calculates conditional correlation.
#' @param counts Single cell expression matrix. Normalizing the data to give every cell the same total expression is preferred.
#' @param conditionon Data frame of variables to be conditioned on when computing gene correlation.
#'   Rows correspond to the rows of counts and xy. At a minimum, it is recommended to include cell type
#'   and tissue ID. Including cell signal strength (total counts) and background (negmean) is also recommended.
#' @param xy Matrix of xy coordinates
#' @param k k for k-nearest neighbor network building
#' @param radius Radius for neighbor network building
#' @param tissue Used for neighbor network building. Neighbors will only be considered for cell with the same tissue value.
#' @param roundcortozero Correlation values with absolute values below this threshold will be stored as 0
#'  to allow for a sparse matrix. If set to NULL, then no rounding to 0 will happen, and a dense matrix will be returned.
#' @param max_cells If there are more than this many cells, certain steps will use a random subset of this size.
#'   Output will still be for all cells.
#' @param return_envmat logical, for whether to return the environment matrix. Default FALSE.
#'   Note the environment matrix is only computed over the subset.
#' @return A list with 3 elements:
#' \itemize{
#' \item condcor: The conditional correlation matrix
#' \item env: The matrix of neighborhood expression of all genes (columns) in all cells - or a subset of cells (rows)
#' \item neighbors: a sparse adjacency matrix giving neighbor relationships
#' }
#' @export
calcSpatialCor <- function(counts, conditionon = NULL,
                       neighbors = NULL, xy = NULL, k = NULL, radius = NULL, tissue = NULL,
                       roundcortozero = 0.1, max_cells = 5000, verbose = TRUE) {

  ## checks: -------------------
  # need xy or neighbors:
  if (is.null(xy) && is.null(neighbors)) {
    stop("need to provide either xy coords or a neighbors network")
  }
  # if xy AND neighbors provided, use neighbors and not xy:
  if (!is.null(xy) && !is.null(neighbors)) {
    message("both xy and neighbors were provided. Using the neighbors provided rather than recalculating.")
    xy <- NULL
  }
  # check the integrity of the neighbors object:
  if (!is.null(neighbors)) {
    if (!identical(dim(neighbors), rep(nrow(counts), 2))) {
      stop("neighbors object must be a matrix or sparse matrix with nrow and ncol both equal nrow(counts)")
    }
    propnonzero <- (sqrt(Matrix::nnzero(neighbors)) / (dim(neighbors)[1]))^2
    if (propnonzero > 0.05) {
      warning("The neighbors matrix provided isn't very sparse - this could indicate an error.")
    }

  }
  # if calculating neighbors, then need k or radius, and not both
  if (!is.null(xy)) {
    if (is.null(k) && is.null(radius)) {
      stop("must provide either k or radius for neighbors to be calculated")
    }
    if (!is.null(k) && !is.null(radius)) {
      message("both k and radius were provided. Proceeding with radius; k will have no impact.")
      k <- NULL
    }
  }
  # warn if no conditionon has been provided:
  if (is.null(conditionon)) {
    warning("conditionon was left NULL. It is recommended to provide at least a vector of cell types.")
  }
  # check format of conditionon:
  if (!is.null(conditionon)) {
    if (nrow(conditionon) != nrow(counts)) {
      stop("conditionon and counts must have the same rows")
    }
  }

  ## build neighbors graph:
  if (is.null(neighbors)) {
    if (verbose) {
      print("building nearest neighbors network")
    }
    if (is.null(tissue)) {
      tissue = 1
    }
    if (!is.null(k)) {
      neighbors <- nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = k, subset = tissue)
    }
    if (!is.null(radius)) {
      neighbors <- radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = radius, subset = tissue)
    }
    rownames(neighbors) <- rownames(xy)
    colnames(neighbors) <- rownames(xy)
  }


  ## perform subsetting if needed:
  use <- TRUE
  if (nrow(counts) > max_cells) {
    use <- sample(seq_len(nrow(counts)), max_cells, replace = FALSE)
  }

  ## get environment matrix: -----------------------
  env <- get_neighborhood_expression(counts = counts,
                                     neighbors = neighbors[use, ])


  ## get the neighborhood values of the conditionon matrix: ------------------- 
  conditionon_neighborhood_vals <- dataframe_neighborhood_summary(df = conditionon, 
                                                                  neighbors = neighbors[use, ])

  ## get matrix to be conditioned on: -------------------------
  condmat <- build_conditional_matrix(conditionon_neighborhood_vals)

  ## get conditional cor: ---------------------
  condcor <- get_conditional_correlation(mat = env, condmat = condmat, outputtype = "cor")
  # round low values to zero to save memory:
  if (!is.null(roundcortozero)) {
    condcor[abs(condcor) < roundcortozero] <- 0
    condcor <- as(condcor, "sparseMatrix")
  }

  ## write output: --------------------------
  out <- list(condcor = round(condcor, 4),
              env = env,
              neighbors = neighbors)
  return(out)
}
