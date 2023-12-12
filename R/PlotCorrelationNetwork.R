#' Plot the conditional correlation matrix as a network
#'
#' Make an igraph plot from a (presumably conditional) correlation matrix. 
#' Makes an adjacency graph between genes by thresholding the correlation matrix,
#' then draws an igraph plot, coloring genes by module membership.#' 
#' @param x Correlation matrix
#' @param modules Either the "modules" data frame output by insitucor(), or a vector of module names.
#' @param genes If \code{modules} is given as a vector, then specify the vector of genes corresponding to it. 
#' @param corthresh Connect genes with cor > this value
#' @param show_gene_names Logical
#' @param vertex_size Argument passed to igraph
#' @param ... Arguments passes to igraph
#' @return A plot.igraph plot
#' @importFrom igraph V
#' @importFrom grDevices colors
#' @importFrom igraph plot.igraph
#' @importFrom uwot umap
#' @export
plotCorrelationNetwork <- function(x, modules, genes = NULL, corthresh = 0.2, show_gene_names = FALSE,
                           vertex_size = NULL, ...) {
  
  # format modules as a named vector:
  if (is.data.frame(modules)) {
    #modules <- setNames(modules$module, modules$gene)
    genes <- modules$gene
    modules <- modules$module
  } 

  # set module colors:
  distinctcolors <-c('#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5', 
                     '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
                     '#e377c2','#17becf','#7f7f7f','#E41A1C','#377EB8',
                     '#BEBADA','#FB8072','#80B1D3','#FDB462','#984EA3',
                     '#A65628','#F781BF','#999999','#FF7F00','#BC80BD',
                     '#8DD3C7','#BEBADA','#FB8072','#80B1D3','#B3DE69',
                     '#BC80BD','#66C2A5','#FC8D62','#8DA0CB','#E78AC3',
                     '#A6D854','#FFD92F','#E5C494','#CCEBC5','#FCCDE5',
                     '#8c564b','#4DAF4A',
                     sample(colors()[!grepl("grey", colors())], 200, replace = FALSE))
  distinctcolors <- distinctcolors[seq_len(length(unique(modules)) + 1)]
  names(distinctcolors) <- c("", unique(modules))
  
  # get a good layout:
  xum <- uwot::umap(as.matrix(x[genes, genes]), 
                   spread = 25, 
                   min_dist = 0.1,
                   n_neighbors = max(min(length(genes)-2, 15), 1))
  #plot(xum, pch = 16, col = distinctcolors[modules])
  
  # make an igraph from the adjacency matrix:
  gr0 <- get_coregulation_network(cormat = x[genes, genes],
                                  corthresh = corthresh)
  
  #  color by module:
  igraph::V(gr0)$module <- modules[match(genes, names(igraph::V(gr0)))]
  
  if (is.null(vertex_size)) {
    vertex_size <- 2 * (1 - show_gene_names)
  }
  igraph::plot.igraph(gr0, 
              layout = xum, 
              vertex.label = unlist(list(NA, names(igraph::V(gr0)))[1 + show_gene_names]), 
              vertex.size = vertex_size, #edge.size = 1, 
              vertex.color = distinctcolors[igraph::V(gr0)$module], 
              vertex.label.color = distinctcolors[igraph::V(gr0)$module],
              ...)
}


#' Get coregulation network
#' 
#' Derive a network graph from a correlation matrix
#' @param cormat The correlation matrix
#' @param corthresh A scalar. Correlations below this value will not produce graph edges.
#' @return An igraph object
#' @importFrom igraph graph_from_adjacency_matrix
get_coregulation_network <- function(cormat, corthresh = 0.3) {
  
  cormatnodiag <- as.matrix(cormat)
  diag(cormatnodiag) <- 0
  cormatnodiag <- round(cormatnodiag, 3)
  
  gr = igraph::graph_from_adjacency_matrix(cormatnodiag * (cormatnodiag > corthresh), weighted = TRUE,
                                           mode = "undirected")
  return(gr)
}
