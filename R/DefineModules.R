
#' Given a conditional correlation matrix, derive modules, including gene weights
#' @param condcor Conditional correlation matrix
#' @param env Environment matrix, possibly for a large subset of cells
#' @param min_module_size Modules smaller than this are discarded
#' @param max_module_size Modules bigger than than this are subclustered
#' @param gene_weighting_rule How to define modules' gene weights, absed on gene expression levels.
#'   One of "inverse_sqrt", "inverse", or "identity".
#' @return A data frame giving module name, gene name, and gene weight, for all genes included in a module.
#' @export
defineModules <- function(condcor, env, min_module_size = 3, max_module_size = 20,
                             resolution = 0.02, corthresh = 0.1, min_module_cor = 0.1,
                             gene_weighting_rule = "inverse_sqrt") {

  modules <- get_modules_from_cor(cormat = condcor,
                                  min_module_size = min_module_size,
                                  max_module_size = max_module_size,             #<----------- max_module_size not implemented yet
                                  resolution = resolution, corthresh = corthresh, min_module_cor = min_module_cor)
  weights <- define_gene_weights(modules = modules,
                                 neighborhood_expression = env,
                                 weighting = gene_weighting_rule)
  # re-format as data frame:
  weightsdf <- rbindlist(lapply(weights, function(x) data.frame(gene = names(x), weight = x)), idcol="module")
  
  # return:
  out = list(modules = modules,
             weights = weights,
             weightsdf = weightsdf)                                             #<----------- would be cleaner if only weightsdf was used downstream
  return(out)
}

#' Get coregulation network using soft thresholded
#'
#' Derive a network graph by soft-thresholding the correlation matrix, and computing Topological Overlap Measure (TOM) values
#' @param cormat The correlation matrix
#' @param min_module_size An integer. Won't consider modules smaller than this.
#' @param max_module_size An integer. Not implemented yet.
#' @param resolution Argument to igraph::cluster_leiden. Lower values produce bigger clusters. 
#' @param corthresh Only correlations about this value will go into the adjacency graph fed into leiden clustering
#' @param min_module_cor Only keep modules with average cor above this value.
#' @return An list with two elements. \code{modules}, a list mod module memberships;
#'  and \code{dend}, a dendrogram from hierarchical clustering of the genes
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_leiden
#' @importFrom igraph E
get_modules_from_cor <- function(cormat, min_module_size = 3, max_module_size = 20,
                                 resolution = 0.02, corthresh = 0.1, min_module_cor = 0.1) {
  
  # get adjacency matrix:
  spcor <- cormat
  spcor <- replace(spcor, spcor < corthresh, 0)
  diag(spcor) = 0
  # convert to a graph, with weights = cor^2
  gr <- igraph::graph_from_adjacency_matrix(adjmatrix = spcor^2 * (spcor>0), 
                                            mode = 'undirected', 
                                            weighted = TRUE)
  # leiden cluster:
  leid <- igraph::cluster_leiden(gr, weights = igraph::E(gr)$weight, resolution_parameter = resolution)
  clust <- leid$membership
  names(clust) <- leid$names
  
  # split excessively large clusters - just once:
  for (cid in unique(clust)) {
    genes <- names(clust)[clust == cid]
    if (length(genes) > max_module_size) {
      # subcluster:
      subgraph <- igraph::graph_from_adjacency_matrix(adjmatrix = spcor[genes, genes]^2 * (spcor[genes, genes] > 0), 
                                                      mode = 'undirected', 
                                                      weighted = TRUE)
      subleid <- igraph::cluster_leiden(subgraph, weights = igraph::E(subgraph)$weight, resolution_parameter = resolution * 2)
      subclust <- subleid$membership
      names(subclust) <- subleid$names
      # replace original cluster id:
      clust[genes] = paste0(clust[genes], "subclust", subclust[genes])
    }
  }
  
  # throw out tiny clusters:
  clustersizes <- table(clust)
  clusternames <- names(clustersizes)[clustersizes >= min_module_size]
  
  # throw out clusters with excessively low correlation, and save the rest in a list:
  modules = list()
  for (cid in clusternames) {
    genes = colnames(cormat)[clust == cid]
    meancor = (sum(cormat[genes, genes]) - length(genes)) / (length(genes)^2 - length(genes)) # direct calculation to not break sparse matrix
    
    if (meancor > min_module_cor) {
      newname <- name_module(cormat[genes, genes])
      modules[[newname]] <- genes
    }
  }
  
  modules = modules[order(sapply(modules, length), decreasing = TRUE)]
  return(modules)
}



#' Score cells for module activity, at single cell level and at neighborhood level
#' @param counts Single cell expression matrix
#' @param weights List of module weights, as output by defineModules
#' @param neighbors Adjacency matrix of neighbors
#' @export
scoreModules <- function(counts, weights, neighbors) {

  ## checks ----------------


  ## compute single cell scores ------------------
  scores_sc <- calculate_metagene_scores(weights = weights, mat = counts)


  ## compute environment scores by spatially averaging the single cell scores --------------------
  scores_env <- get_neighborhood_expression(counts = scores_sc,                    
                                            neighbors = neighbors)

  ## write output -----------
  out <- list(scores_env = scores_env,
              scores_sc = scores_sc)
  return(out)
}

#' Define weights for cluster metagenes
#'
#' Given a cluster's gene IDs and the environment expression matrix, define how to weight genes for cluster metagene scores
#' @param modules A list, with each element a vector of gene names
#' @param neighborhood_expression The neighborhood expression matrix. Colnames must contain the genes given in clusts
#' @param weighting One of "inverse_sqrt", "inverse", or "identity"
#' @return A list, where each element is a named vector giving the gene weights for a cluster.
define_gene_weights <- function(modules, neighborhood_expression, weighting = "inverse_sqrt") {
  if (!is.element(weighting, c("inverse_sqrt", "inverse", "identity"))) {
    stop("weighting argument must be one of \'inverse_sqrt\', \'inverse\', or \'identity\'")
  }
  meanneighborhoodexpr <- Matrix::colMeans(neighborhood_expression)

  wts <- list()
  for (i in 1:length(modules)) {
    if (weighting == "inverse_sqrt") {
      wts[[i]] <- meanneighborhoodexpr[modules[[i]]]^-0.5
    }
    if (weighting == "inverse") {
      wts[[i]] <- meanneighborhoodexpr[modules[[i]]]^-1
    }
    if (weighting == "identity") {
      wts[[i]] <- rep(1, length(modules[[i]]))
      names(wts[[i]]) <- modules[[i]]
    }
    wts[[i]] <- wts[[i]] / sum(wts[[i]])
  }
  names(wts) <- names(modules)
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
  scores <- as.matrix(scores)
  return(scores)
}


#' name module based on gene PCs
#'
#' name a module based on its top 2 genes
#' @param mat The conditional correlation matrix for the selected genes
#' @return a name
#' @importFrom Matrix colMeans
name_module <- function(mat) {
  meancors <- Matrix::colMeans(mat)
  n <- ncol(mat)
  top3 <- colnames(mat)[order(meancors, decreasing = TRUE)[1:min(n, 3)]]
  if (n <= 3) {
    name <- paste0(c(top3, n), collapse = "_")
  } else {
    name <- paste0(c(top3[1:2], n), collapse = "_")
  }
  return(name)
}

