



#### Reqs for sparc:
SPARC calls the complete workflow for spatial correlation analysis. 

##### Inputs:
- counts Single cell expression matrix. Normalizing the data to give every cell the same total expression is preferred.
- conditionon Data frame of variables to be conditioned on when computing gene correlation. Rows correspond to the rows of counts and xy. At a minimum, it is recommended to include cell type and tissue ID. Including cell signal strength (total counts) and background (negmean) is also recommended.
- celltype Vector of cell type assignments
- neighbors Neighbor relationships, stored as a sparse matrix
- xy Matrix of xy coordinates.
- k k for k-nearest neighbor network building
- radius Radius for neighbor network building
- tissue Used for neighbor network building. Neighbors will only be considered for cell with the same tissue value.
- mincor Correlation values below this threshold will be stored as 0 to allow for a sparse matrix
- resolution Argument to igraph::cluster_leiden. Lower values produce bigger clusters. 
- corthresh Only correlations about this value will go into the adjacency graph fed into leiden clustering
- min_module_cor Only keep modules with average cor above this value.
- max_cells If there are more than this many cells, certain steps will use a random subset of this size.


##### Outputs:
A list, with the following elements:
\enumerate{
  \item condcor: A sparse matrix holding genes' conditional correlations.
  \item modules: A data frame detailing gene module membership and weights
  \item scores_env: A matrix of cell * module environment scores
  \item scores_sc: A matrix of cell * module single cell scores
  \item attributionmats: A list of matrices holding attribution scores for each cell type * gene in each module.
  \item celltypeinvolvement: A matrix giving the maximum attribution score for each cell type in each module.
    }



#### Reqs for calcSpatialCor:
All the steps to compute the conditional spatial correlation matrix.

##### Inputs:
- counts Single cell expression matrix. Normalizing the data to give every cell the same total expression is preferred.
- conditionon Data frame of variables to be conditioned on when computing gene correlation. Rows correspond to the rows of counts and xy. At a minimum, it is recommended to include cell type and tissue ID. Including cell signal strength (total counts) and background (negmean) is also recommended.
- neighbors Neighbor relationships, stored as a sparse matrix
- xy Matrix of xy coordinates
- k k for k-nearest neighbor network building
- radius Radius for neighbor network building
- tissue Used for neighbor network building. Neighbors will only be considered for cell with the same tissue value.
- roundcortozero Correlation values with absolute values below this threshold will be stored as 0
  to allow for a sparse matrix. If set to NULL, then no rounding to 0 will happen, and a dense matrix will be returned.
- max_cells If there are more than this many cells, certain steps will use a random subset of this size. Output will still be for all cells.
- verbose Whether to print progress

##### Outputs:
A list, with the following elements:
\enumerate{
\item condcor: The conditional correlation matrix
\item env: The matrix of neighborhood expression of all genes (columns) in all cells - or a subset of cells (rows)
\item neighbors: a sparse adjacency matrix giving neighbor relationships
  }




#### Reqs for defineModules:
Define gene modules given the spatial correlation matrix.

##### Inputs:
- condcor Conditional correlation matrix
- env Environment matrix, possibly for a large subset of cells
- min_module_size Modules smaller than this are discarded
- max_module_size Modules bigger than than this are subclustered
- gene_weighting_rule How to define modules' gene weights, absed on gene expression levels. One of "inverse_sqrt", "inverse", or "identity".
- resolution Resolution parameter for leiden clustering
- corthresh Correlations with absolute value below this will be rounded to zero to save memory
- min_module_cor Modules must have mean correlation of at least this much to be reported


##### Outputs:
A list, with the following elements:
\enumerate{
\item modules: A list giving the gene names in each module
\item weights: A list giving the gene weights in each module
\item weightsdf: A data frame summary of genes' module membership and weights
  }
  
  
  
#### Reqs for scoreModules:
Score gene modules over single cells and their environments.

##### Inputs:
- counts Single cell expression matrix
- weights List of module weights, as output by defineModules
- neighbors Adjacency matrix of neighbors


##### Outputs:
A list, with the following elements:
\enumerate{
\item scores_env: A matrix of module scores for cell environments
\item scores_sc: A matrix of module scores for single cell profiles
  }
  
  
#### Reqs for cellTypeAttribution:
Score cell types' involvement in modules
##### Inputs:
- modulescores Matrix of each cell's score for each module.
- weights List of module gene weights.
- counts Single cell expression matrix, cells in rows, genes in columns. The same data used in earlier steps.
- celltype Vector of cell types
- neighbors Neighbor relationships, stored as a sparse matrix
- nsub Subsample size to use. This function only extracts correlations, so 5000 cells is ample
- verbose Whether to print progress


##### Outputs:

A list, with the following elements:
\enumerate{
\item "involvescores", a matrix of module * cell type scores, from 0-1.
\item "attributionmats", a list of matrices giving the involvement of each cell type in each gene in each cluster.
  }
  
  
  
#### Reqs for plotCorrelationNetwork:
Make an igraph plot of the network structure of the conditional correlation matrix:
##### Inputs:
- x Correlation matrix
- modules Either the "modules" data frame output by sparc(), or a vector of module names.
- genes If modules is given as a vector, then specify the vector of genes corresponding to it. 
- corthresh Connect genes with cor > this value
- show_gene_names Logical
- vertex_size Argument passed to igraph


##### Outputs:

Draws an igraph network plot
  
  