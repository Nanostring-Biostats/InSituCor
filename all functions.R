#library(WGCNA)

#' Normalize raw counts for use in method
#'
#' Background-subtracts and scales raw counts
#' @param raw Matrix of raw counts. Cells in rows.
#' @param neg Vector giving each cell's mean negprobe expression
#normalize_for_neighborhood_analysis(raw, neg) {
#
#}

#' Calculate neighborhood expression profiles
#'
#' Calculates the total expression in the neighborhood around each cell
#' @param counts Matrix of counts, with cells in rows and genes in columns. Can be normalized, but must be linear-scale. (I.e. cannot be log-scale or Pearson residuals)
#' @param xy Matrix of xy locations, with rows matching counts
#' @param k Number of nearest neighbors to use
#' @param neighborsmatrix Matrix giving the row numbers of each cell's k nearest neighbors.
#'  If NULL, then xy and k will be used to derive this matrix.
#'  Inputting this matrix will save computation time.
#' @param neg Vector giving each cell's mean negprobe expression. If provided then neighborhood expression will be background-subtracted.
#' @param scaling Vector giving scaling factors for each cell. If provided then each cell's neighborhood expression will be scaled by the neighborhood total of this factor.
#' @return A matrix giving the neighborhood expression around each cell.
get_neighborhood_expression <- function(counts, xy=NULL, k=50, neighborsmatrix=NULL, verbose = FALSE, neg = NULL, scaling = NULL) {
  if (is.null(xy) & is.null(neighborsmatrix)) {
    stop("need either xy or neighborsmatrix")
  }
  # get k nearest neighbors in xy space
  if (is.null(neighborsmatrix)) {
    if (verbose) {
      print("building nearest neighbors matrix")
    }
    neighborsmatrix <- spatstat.geom::nnwhich(X = xy[, 1], Y = xy[, 2], k = 1:k)
    rownames(neighborsmatrix) <- rownames(xy)
  }
  # get total neighborhood expression around each cell:
  totalneighborhoodexpr = as.matrix(counts * NA)
  if (verbose) {
    print("calculating total expression over neighbors")
  }
  for (i in 1:ncol(counts)) {
    if (verbose & (i %% 50 == 0)) {print(i)}
    totalneighborhoodexpr[, i] <- rowSums(matrix(counts[, i][as.vector(neighborsmatrix)], nrow = nrow(neighborsmatrix)))
  }
  # background subtract and normalize:
  if (!is.null(neg)) {
    totalneighborhoodneg <- rowSums(matrix(annot$negmean[as.vector(neighborsmatrix)], nrow = nrow(neighborsmatrix)))
    totalneighborhoodexpr <- pmax(sweep(totalneighborhoodexpr, 1, totalneighborhoodneg, "-"), 0)
  }
  if (!is.null(scaling)) {
    totalneighborhoodexpr <- sweep(totalneighborhoodexpr, 1, scaling, "/")
  }
  return(totalneighborhoodexpr)
}


#' Calculate neighborhood cell type makeup
#'
#' Calculates the total expression in the neighborhood around each cell
#' @param celltypes Vector of cell types
#' @param xy Matrix of xy locations, with rows aligned to celltypes
#' @param k Number of nearest neighbors to use
#' @return A matrix giving the cell type abundance in each cell's neighborhood
get_neighborhood_celltypes <- function(celltypes, xy, k=50) {
  pointprocess <- spatstat.geom::ppp(xy[, 1], xy[, 2],
                                     range(xy[, 1]), range(xy[, 2]),
                                     marks = as.factor(celltypes),
                                     unitname = c("mm", "mm"))
  celltypeneighborhoods <- spatstat.core::marktable(pointprocess, N = k, exclude=FALSE)
  return(celltypeneighborhoods)
}

#' Combine all variables to be conditioned on
#' @param variables A list or data frame containing all the variables (could be vectors or matrices) that you wish to condition on
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


##' evaluate a matrix for whether any variables will pose problems for inversion:
#flag_problematic_variables <- function(mat) {
#
#}

#' get conditional covariance of one matrix given another:
#' @param mat Matrix for which conditional covariance is desired
#' @param condmat Matrix to be conditioned on
#' @param returncor Logical, for whether to convert to correlation
#' @return A conditional covariance (or correlation) matrix
get_conditional_covariance <- function(mat = NULL, condmat = NULL, returncor = FALSE) {
  cov_aa <- cov(mat)
  cov_bb <- cov(condmat)
  cov_ab <- cov(mat, condmat)

  condcov <- cov_aa - cov_ab %*% solve(cov_bb) %*% t(cov_ab)

  if (returncor) {
    condcov <- cov2cor(condcov)
  }
  return(condcov)
}

#' convert covariance matrix to correlation:
cov2cor <- function(covmat) {
  vars <- diag(covmat)
  cormat <- diag(vars^-0.5) %*% covmat %*% diag(vars^-0.5)
  dimnames(cormat) <- dimnames(covmat)
  return(cormat)
}


#' convenience function to wrap up the below 2 functions
#
#  FUNCTION HERE () {
#
#
#
#
#  }


#' Get coregulation network
#'
#' Derive a network graph from a correlation matrix
#' @param cormat The correlation matrix
#' @param corthresh A scalar. Correlations below this value will not produce graph edges.
#' @return An igraph object
#' @importFrom igraph graph_from_adjacency_matrix
get_coregulation_network <- function(cormat, corthresh = 0.3) {

  cormatnodiag = cormat
  diag(cormatnodiag) = 0

  gr = igraph::graph_from_adjacency_matrix(cormatnodiag * (cormatnodiag > corthresh), weighted = TRUE,
                                           mode = "undirected")
  return(gr)
}


#' Cluster coregulation network
#'
#' Cluster a network built from a correlation matrix, using leiden clustering
#' @param gr An igraph object
#' @param resolution_parameter The resolution parameter for cluster_leiden()
#' @param ... Other parameters passed to cluster_leiden
#' @return
#' @importFrom igraph cluster_leiden
cluster_coreg_net <- function(gr, resolution_parameter = 0.2, min_cluster_size = 2, ...) {

  # get leiden clusters:
  cl <- igraph::cluster_leiden(gr, resolution_parameter = .2, objective_function = "CPM")

  # save list of cluster memberships:
  modules <- list()
  for (id in names(which(table(cl$membership) >= min_cluster_size))) {
    modules[[id]] <- cl$names[cl$membership == id]
  }
  modules <- modules[order(sapply(modules, length), decreasing = TRUE)]
  names(modules) <- paste0("m", seq_len(length(modules)))
  return(modules)
}


#' name module based on gene PCs
#'
#' name a module based on its top 2 genes
#' @param mat The conditional correlation matrix for the selected genes
#' @return a name
name_module <- function(mat) {
  meancors <- colMeans(mat)
  n <- ncol(mat)
  top3 <- colnames(mat)[order(meancors, decreasing = TRUE)[1:min(n, 3)]]
  if (n <= 3) {
    name <- paste0(c(top3, n), collapse = "_")
  } else {
    name <- paste0(c(top3[1:2], n), collapse = "_")
  }
  return(name)
}




#' Get coregulation network using soft thresholded
#'
#' Derive a network graph by soft-thresholding the correlation matrix, and computing Topological Overlap Measure (TOM) values
#' @param cormat The correlation matrix
#' @param min_module_size An integer. Won't consider modules smaller than this.
#' @param max_module_size An integer. Not implemented yet.
#' @param beta A scalar. Soft thresholding is (0.5 + 0.5*cor)^beta
#' @param deepSplit An integer in 0:4, used by dynamicTreeCut. Smaller produces more genes in clusters and bigger clusters.
#' @param
#' @return An igraph object
#' @importFrom igraph graph_from_adjacency_matrix
get_modules_w_WGCNA <- function(cormat = condcor, min_module_size = 3, max_module_size = 20,
                                beta = 12, deepSplit = 4, mincor = 0.2) {

  # get TOM distance:
  softthresh <- (0.5 + 0.5 * cormat)^beta
  #disssoft=1-softthresh
  #simTOM=TOMsimilarity(softthresh)
  dissTOM=TOMdist(softthresh)
  dimnames(dissTOM) = dimnames(softthresh)
  collectGarbage()

  hierTOM = hclust(as.dist(dissTOM),method="average")

  # get a dynamic tree cut:
  cuts = dynamicTreeCut::cutreeDynamic(dendro = hierTOM, distM = dissTOM,
                                       minClusterSize = min_module_size,
                                       method = "hybrid", pamStage = FALSE, deepSplit = deepSplit)

  # throw out cuts with excessively low correlation, and save the rest in a list:
  unclusteredvalue = 0
  modules = list()
  for (cid in setdiff(unique(cuts), unclusteredvalue)) {
    genes = colnames(cormat)[cuts == cid]
    meancor = mean(cormat[genes, genes][upper.tri(cormat[genes, genes])])
    if (meancor > mincor) {
      newname <- name_module(cormat[genes, genes])
      modules[[newname]] <- genes
    }
  }
  modules = modules[order(sapply(modules, length), decreasing = TRUE)]
  out = list(modules = modules, dend = hierTOM)
  return(out)
}



#' Define layout for a coregulation network
#'
#' Default plot for a coregulation network. Uses the igraph layout "
#' @param gr igraph object
#' @param showgenes Vector of gene names to show. If NULL, then all connected genes will be shown.
#layout_coreg_net <- function() {
#
#}


#' Plot a coregulation network
#'
#' Default plot for a coregulation network. Uses the igraph layout "
#' @param gr igraph object
#' @param genes Vector of gene names for genes to show. If NULL, then all connected genes will be shown.
#' @param modules Named vector of cluster memberships
#' @param show_gene_names Logical
#' @param genes_to_show If given a vector of gene names, will only show these gene names
#' @param ... Arguments passed to plot.igraph
plot_coreg_net <- function(gr, genes = NULL, modules = NULL, show_gene_names = FALSE, genes_to_show = NULL, ...) {
  distinctcolors <- c('#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5',   #'#c7c7c7',
                               '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b',
                               '#e377c2','#17becf','#7f7f7f',
                               '#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462',
                               '#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F',
                               '#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5',
                               '#D9D9D9','#BC80BD','#CCEBC5','#FFED6F','#66C2A5','#FC8D62','#8DA0CB','#E78AC3',
                               '#A6D854','#FFD92F','#E5C494','#B3B3B3','#E41A1C','#377EB8','#4DAF4A','#984EA3',
                               '#FF7F00','#FFFF33','#A65628','#F781BF','#999999',
                               sample(colors()[!grepl("grey", colors())], 200, replace = FALSE))
  distinctcolors <- distinctcolors[seq_len(length(modules) + 1)]
  names(distinctcolors) <- c("", names(modules))

  # subset graph:
  if (is.null(genes)) {
    # get all connected genes:
    genes <- names(which(degree(gr) > 0))
  }
  gr0 <- delete_vertices(gr, setdiff(names(V(gr)), genes))

  # label vertices with modules:
  modulenames <- rep("", length(V(gr0)))
  names(modulenames) <- names(V(gr0))
  if (!is.null(modules)) {
    for (i in seq_len(length(modules))) {
      modulenames[modules[[i]]] <- names(modules)[i]
    }
  }
  #  color by module:
  V(gr0)$module = modulenames[names(V(gr0))]

  # keep only the desired gene labels:  (this didn't work)
  #if (is.null(genes_to_show)) {
  #  genes_to_show <- names(V(gr0))
  #}
  #names_to_show <- rep("", length(V(gr0)))
  #names_to_show[match(genes_to_show, names(V(gr0)))] <- genes_to_show
  #names(V(gr0)) <- names_to_show

  # layout graph and plot:
  layout <- layout_with_lgl(gr0)
  igraph::plot.igraph(gr0, coords = layout,
                      vertex.label = unlist(list(NA, names(V(gr0)))[1 + show_gene_names]),
                      vertex.size = 2 * (1 - show_gene_names), edge.size = 1,
                      #vertex.color = distinctcolors[modulenames[names(V(gr0))]],
                      vertex.color = distinctcolors[V(gr0)$module],
                      vertex.label.color = distinctcolors[V(gr0)$module],
                      ...)
}


#' Define weights for cluster metagenes
#'
#' Given a cluster's gene IDs and the environment expression matrix, define how to weight genes for cluster metagene scores
#' @param clusts A list, with each element a vector of gene names
#' @param neighborhood_expression The neighborhood expression matrix. Colnames must contain the genes given in clusts
#' @param weighting One of "inverse_sqrt", "inverse", or "identity"
#' @return A list, where each element is a named vector giving the gene weights for a cluster.
derive_gene_weights <- function(clusts, neighborhood_expression, weighting = "inverse_sqrt") {
  if (!is.element(weighting, c("inverse_sqrt", "inverse", "identity"))) {
    stop("weighting argument must be one of \'inverse_sqrt\', \'inverse\', or \'identity\'")
  }
  meanneighborhoodexpr <- Matrix::colMeans(neighborhood_expression)

  wts <- list()
  for (i in 1:length(clusts)) {
    if (weighting == "inverse_sqrt") {
      wts[[i]] <- meanneighborhoodexpr[clusts[[i]]]^-0.5
    }
    if (weighting == "inverse") {
      wts[[i]] <- meanneighborhoodexpr[clusts[[i]]]^-1
    }
    if (weighting == "identity") {
      wts[[i]] <- rep(1, length(clusts[[i]]))
      names(wts[[i]]) <- clusts[[i]]
    }
    wts[[i]] <- wts[[i]] / sum(wts[[i]])
  }
  names(wts) <- names(clusts)
  return(wts)
}


#' Calculate metagene scores
#'
#' @param weights List of named vectors giving gene weights, such as is output by derive_gene_weights()
#' @param mat Expression matrix, cells in rows and genes in columns. Could be either environment or single cell expression.
calculate_metagene_scores <- function(weights, mat) {
  scores <- c()
  for (i in 1:length(weights)) {
    scores <- cbind(scores, mat[, names(weights[[i]])] %*% weights[[i]])
  }
  colnames(scores) <- names(weights)
  return(scores)
}



#' score modules for celltype involvement
#'
#' Each module * cell type is given a score, definined by the correlation between environment score
#'  and cell type specific environment gene expression cell type's most
#' impressive contribution to a gene's expression in module hotspots.
#' @param modules List of module memberships.
#' @param scores Matrix of module scores.
#' @param counts Single cell expression matrix, cells in rows, genes in columns
#' @param clusts Vector of cell types
#' @param neighborsmatrix Matrix assigning cells to neighbors
#' @param neg Vector of negprobes
#' @param scaling Vector of scaling factors
#' @param nsub Subsample size to use. This function only extracts correlations, so 5000 cells is ample
#' @return A list with two elements: "involvescores", a matrix of module * cell type scores, from 0-1.
#'         "attributionmats", a list of matrices giving the involvement of each cell type in each gene in each cluster.
score_module_celltype_involvement <- function(modules, scores, counts, clusts, neighborsmatrix,
                                              neg = NULL, scaling = NULL, nsub = 5000, verbose = FALSE) {
  celltypes <- unique(clusts)
  involvescores <- matrix(0, length(modules), length(celltypes),
                          dimnames = list(names(modules), celltypes))
  attributionmats <- list()

  for (name in names(modules)) {
    if (verbose) {
      message(name)
    }
    # identify hotspots:
    attributionmats[[name]] <- get_single_module_celltype_involvement(
      score_env = scores[, name],
      counts = counts,
      genes = modules[[name]],
      clusts = clusts,
      neighborsmatrix = neighborsmatrix,
      neg = neg,
      scaling = scaling,
      nsub = nsub)
    involvescores[name, rownames(attributionmats[[name]])] <- apply(attributionmats[[name]], 1, max)
  }
  out <- list(involvescores = involvescores, attributionmats = attributionmats)
  return(out)
}




#' score correlation between environment score and per-cell type expression in neighborhood
#'
#' NOTE: this one works great and is fast enough. Still needs a wrapper in the vein of "score_module_celltype_involvement()".
#' @param score Vector of environment scores
#' @param counts Counts matrix, cells in row, columns holding the genes given in the "genes" argument
#' @param genes Vector of gene names
#' @param clusts Vector of cell type assignments
#' @param neighborsmatrix Matrix assigning cells to neighbors
#' @param neg Vector of negprobes
#' @param scaling Vector of scaling factors
#' @param nsub Subsample size to use. This function only extracts correlations, so 5000 cells is ample
get_single_module_celltype_involvement <- function(score_env, counts, genes, clusts, neighborsmatrix,
                                                   neg = NULL, scaling = NULL, nsub = 5000) {

  ## subsampling details:
  if (is.null(nsub)) {
    use <- TRUE
  } else {
    use <- sample(1:length(score_env), min(nsub, length(score_env)))
  }
  ## result to return: attribution stat for each clust * gene:
  attributionmat <- matrix(NA, nrow = length(unique(clusts)), ncol = length(genes),
                           dimnames = list(unique(clusts), genes))

  ## build cell type specific env matrices for these genes:
  clust_envmats <- list()
  for (cell in unique(clusts)) {
    clust_envmats[[cell]] <- matrix(0, length(score_env[use]), length(genes),
                                    dimnames = list(names(score_env)[use], genes))
  }

  for (cell in unique(clusts)) {
    #print(cell)

    # record environment clust-specific neg for subtraction:
    if (!is.null(neg)) {
      totalneighborhoodneg <- rowSums(matrix(
        (annot$negmean * (clusts == cell))[as.vector(neighborsmatrix[use, ])], nrow = nrow(neighborsmatrix[use, ])))
    }

    # get clust-specific environment expression
    for (gene in genes) {
      clust_envmats[[cell]][, gene] <-
        rowSums(matrix((counts[, gene] * (clusts == cell))[as.vector(neighborsmatrix[use, ])], nrow = nrow(neighborsmatrix[use, ])))
      if (!is.null(neg)) {
        clust_envmats[[cell]][, gene] <- pmax(clust_envmats[[cell]][, gene] - totalneighborhoodneg, 0)
      }
    }

    # apply same scaling factors used by "env":
    if (!is.null(scaling)) {
      clust_envmats[[cell]] <- sweep(clust_envmats[[cell]], 1, scaling[use], "/")
    }

    # record correlation w score:
    attributionmat[cell, genes] <- cor(score_env[use], clust_envmats[[cell]], use = "pairwise.complete")
  }
  attributionmat <- replace(attributionmat, is.na(attributionmat), 0)
  ## summarize correlation between matrices and
  #out = list(attributionmat = attributionmat, clust_envmats = clust_envmats)
  #return(out)
  return(attributionmat)
}




#### older versions, no longer used -----------------------------------

if (FALSE) {

  #' score modules for celltype involvement
  #'
  #' Each module * cell type is given a score, definined by that cell type's most
  #' impressive contribution to a gene's expression in module hotspots.
  #' @param modules List of module memberships.
  #' @param scores Matrix of module scores.
  #' @param counts Single cell expression matrix, cells in rows, genes in columns
  #' @param clusts Vector of cell types
  #' @param hotspotquant Scores about this quantile will be considered hotspots. (NOTE: this could be refined; it results in an arbitrary hotspot definition)
  #' @return A list with two elements: "involvescores", a matrix of module * cell type scores, from 0-1.
  #'         "attributionmats", a list of matrices giving the involvement of each cell type in each gene in each cluster.
  score_module_celltype_involvement_old <- function(modules, scores, counts, clusts, hotspotquant = 0.99) {
    celltypes <- unique(clusts)
    involvescores <- matrix(0, length(modules), length(celltypes),
                            dimnames = list(names(modules), celltypes))
    attributionmats <- list()

    for (name in names(modules)) {
      genes <- modules[[name]]

      # identify hotspots:
      inscorebin <- scores[, name] > quantile(scores[, name], hotspotquant)
      temp <- by(as.matrix(counts[inscorebin, genes, drop = FALSE]), clusts[inscorebin], colMeans)
      temp <- sapply(temp, cbind)
      rownames(temp) <- genes
      temp <- sweep(temp, 1, rowSums(temp), "/")
      attributionmats[[name]] <- temp
      tempcellscores <- apply(temp, 2, max)
      involvescores[name, names(tempcellscores)] <- tempcellscores
    }
    out <- list(involvescores = involvescores, attributionmats = attributionmats)
    return(out)
  }



  #' Attribute module genes to cell types
  #'
  #' For each module, get the amount of expression each gene is getting from each cell type within module "hotspots"
  #' @param modules List of module memberships
  #' @param scores Matrix of model scores
  #' @param celltype Vector of cell types, aligned to rows of scores
  #' @param counts Matrix of single cell gene expression, aligned to rows of scores
  #' @param quant Scores above this threshold will be considered "hotspots"
  #' @return A list of matrices, with each matrix giving the cell type * gene expression levels in hotspots
  get_module_celltype_involvement_old <- function(modules, scores, celltype, counts, quant = 0.95) {
    attributionmats <- list()
    for (cl in names(modules)) {
      genes <- modules[[cl]]
      inscorebin <- scores[, cl] > quantile(scores[, cl], quant)
      # get total expression per cell type:
      #sumbycelltype <- matrix(0, length(genes), length(unique(celltype)),
      #                        dimnames = list(genes, unique(celltype)))
      #for (cell in unique(celltype)) {
      #  sumbycelltype[, cell] <- colSums(counts[inscorebin & (celltype == cell), genes, drop = FALSE])
      #}
      #attributionmats[[cl]] <- sumbycelltypes
      temp <- by(as.matrix(counts[inscorebin, genes, drop = FALSE]), clusts[inscorebin], colMeans)
      temp <- sapply(temp, cbind)
      rownames(temp) <- genes
      attributionmats[[cl]] <- temp
    }
    return(attributionmats)
  }

}
