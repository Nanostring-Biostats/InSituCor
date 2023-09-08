
#' score modules for celltype involvement
#'
#' Each module * cell type is given a score, defined by the correlation between environment score
#'  and cell type specific environment gene expression cell type's most
#' impressive contribution to a gene's expression in module hotspots.
#' @param modulescores Matrix of each cell's score for each module.
#' @param weights List of module gene weights.
#' @param counts Single cell expression matrix, cells in rows, genes in columns. The same data used in earlier steps.
#' @param celltype Vector of cell types
#' @param neighbors Neighbor relationships, stored as a ___________
#' @param nsub Subsample size to use. This function only extracts correlations, so 5000 cells is ample
#' @return A list with two elements: "involvescores", a matrix of module * cell type scores, from 0-1.
#'         "attributionmats", a list of matrices giving the involvement of each cell type in each gene in each cluster.
#' @export
cellTypeAttribution <- function(modulescores, weights, counts, celltype, neighbors,
                                 nsub = 5000, verbose = FALSE) {
  # initialize:
  if(is.factor(celltype)) {
    celltype <- as.character(celltype)
  }
  celltypes <- unique(celltype)
  involvescores <- matrix(0, length(weights), length(celltypes),
                          dimnames = list(names(weights), celltypes))
  attributionmats <- list()

  ## perform subsampling:
  if (is.null(nsub)) {
    use <- TRUE
  } else {
    use <- sample(1:nrow(counts), min(nsub, nrow(counts)))
  }
  modulescores <- modulescores[use, ]
  neighbors <- neighbors[use, ]

  # compute cell type attribution stat for each module, using subsampled data:
  for (name in names(weights)) {
    if (verbose) {
      message(name)
    }
    # score cell type x gene attribution for a single module
    attributionmats[[name]] <- get_single_module_celltype_involvement(
      scorevector = modulescores[, name],
      genes = names(weights[[name]]),
      counts = counts[, names(weights[[name]])],
      celltype = celltype,
      neighbors = neighbors)
    # add columns for rare cell types that didn't get considered:
    missingcelltypes <- setdiff(celltypes, rownames(attributionmats[[name]]))
    if (length(missingcelltypes) > 0) {
      missingattribmat <- matrix(0, length(missingcelltypes), ncol(attributionmats[[name]]),
                                 dimnames = list(missingcelltypes, colnames(attributionmats[[name]])))
      attributionmats[[name]] <- rbind(attributionmats[[name]], missingattribmat)
    }
    attributionmats[[name]] <- attributionmats[[name]][celltypes, ]

    # summarize each cell type's result for the module:
    involvescores[name, rownames(attributionmats[[name]])] <- apply(attributionmats[[name]], 1, max)
  }
  out <- list(involvescores = involvescores, attributionmats = attributionmats)
  return(out)
}




#' score correlation between environment score and per-cell type expression in neighborhood
#'
#' For all cell types x genes in the module, reports the correlation between module score and environment expression of the 
#' gene coming from the cell type.
#' @param scorevector Vector giving each cell's environment score for the module. Possibly a subset.
#' @param genes Vector of module gene names
#' @param counts Counts matrix, cells in row, columns holding the genes given in the "genes" argument. Not subsetted. 
#' @param celltype Vector of cell type assignments. Not subsetted.
#' @param neighbors Neighbor relationships, stored as a sparse matrix. Possibly rows are just a subset, but columns are present for all cells.
get_single_module_celltype_involvement <- function(scorevector, genes, counts, celltype, neighbors) {


  ## result to return: attribution stat for each clust * gene:
  attributionmat <- matrix(NA, 
                           nrow = length(unique(celltype)), 
                           ncol = length(genes),
                           dimnames = list(unique(celltype), genes))

  ## for each cell type, calc a matrix of environment expression attributable to the cell type:
  for (cell in unique(celltype)) {
    if (sum(celltype == cell) > 10) {
      
      # get neighborhood expression attributable to the cell type:
      neighborhood_expression_in_celltype <- get_neighborhood_expression(
        counts = sweep(counts[, genes], 1, 1 * (celltype == cell), "*"), 
        neighbors = neighbors)
      
      # save correlation between module score and cell type-specific neighborhood expression:
      attributionmat[cell, ] <- suppressWarnings(
        cor(scorevector, 
            neighborhood_expression_in_celltype[, genes, drop = FALSE],
            use = "pairwise.complete"))
      attributionmat[cell, ] <- replace(attributionmat[cell, ], is.na(attributionmat[cell, ]), 0)
    }
  }
  
  attributionmat <- replace(attributionmat, is.na(attributionmat), 0)
  return(attributionmat)
}
