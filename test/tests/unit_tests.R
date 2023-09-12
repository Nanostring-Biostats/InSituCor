library(testthat)
library(InSituCor)
rm(list = ls())
#temp <- readRDS("data/mini CosMx kidney.RDS")
#annot <- temp$annot
#counts <- temp$counts
#celltype <- temp$celltype
#xy <- temp$xy
#rm(temp)
data(cosmx_kidney)
annot <- cosmx_kidney$annot
rownames(annot) <- annot$cell_ID
counts <- cosmx_kidney$counts
celltype <- as.factor(cosmx_kidney$annot$celltype)
xy <- cosmx_kidney$xy
neighbors <- InSituCor:::radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = 0.05, subset=1)


head(rownames(counts))
head(rownames(annot))

# for running functions line by line:
if (FALSE) {
  conditionon = cbind(annot[, c("fov", "totalcounts")], celltype)
  neighbors = NULL
  k = 50
  radius = 0.05
  tissue = NULL
  mincor = 0.1
  max_cells = 1e5
  verbose = TRUE
}


#### test big wrapper -----------------------------------

res <- insitucor(counts = counts,
             conditionon = annot[, c("fov", "totalcounts", "celltype")],
             celltype = annot$celltype,
             neighbors = NULL, xy = xy, k = NULL, radius = 0.05, tissue = annot$fov, # args for neighbor definition
             min_module_size = 2, max_module_size = 8,                 # args for module definition
             gene_weighting_rule = "inverse_sqrt",   # more args for module definition
             roundcortozero = 0.1, max_cells = 1e5,                               # args for controlling memory and compute
             attribution_subset_size = 1000,                                      # args for cell type attribution scoring
             verbose = TRUE)

test_that("wrapper returns the expected results", {
  expect_identical(colnames(res$modules), c("module", "gene", "weight"))
  expect_true(nrow(res$scores_env) == nrow(counts))
  expect_true(nrow(res$scores_sc) == nrow(counts))
  expect_true(is(res$condcor, 'sparseMatrix'))
  expect_true(is.matrix(res$celltypeinvolvement))
  expect_true(is.matrix(res$attributionmats[[1]]))
})


#### test neighbor definition: -----------------------

neighbors <- InSituCor:::radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = 0.05, subset=1)
neighbors <- InSituCor:::nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = 50, subset=1)
# (should be sparse matrices with dim = nrow(xy).)


#### test neighbor summing math: --------------------

test_that("neighbor math has the right logic", {
  
  expect_equal(neighbor_tabulate(annot$celltype, neighbors)[1:5, 1:3],
               neighbor_tabulate(annot$celltype, neighbors[1:5,])[, 1:3])
  
  expect_equal(neighbor_sum(annot$totalcounts, neighbors)[1:5],
               neighbor_sum(annot$totalcounts, neighbors[1:5, ]))
  
  expect_equal(neighbor_colSums(counts, neighbors)[1:5, ],
               neighbor_colSums(counts, neighbors[1:5, ]))
  
  expect_equal(neighbor_meantabulate(annot$celltype, neighbors)[1:5, 1:3],
               neighbor_meantabulate(annot$celltype, neighbors[1:5,])[, 1:3])
  
  expect_equal(neighbor_mean(annot$totalcounts, neighbors)[1:5],
               neighbor_mean(annot$totalcounts, neighbors[1:5, ]))
  
  expect_equal(neighbor_colMeans(counts, neighbors)[1:5, ],
               neighbor_colMeans(counts, neighbors[1:5, ]))
 
  expect_equal(neighbor_colMeans(counts, neighbors)[1:5, 1],
               neighbor_mean(counts[, 1], neighbors)[1:5])
 
  expect_equal(neighbor_colSums(counts, neighbors)[1:5, 1],
               neighbor_sum(counts[, 1], neighbors)[1:5])
  
  expect_equal(neighbor_colSums(counts, neighbors)[1:5, 1],
               neighbor_colSums(as.matrix(counts), neighbors)[1:5, 1])
})


#### test functions calling neighborhood summaries: ----------------
tmp <- dataframe_neighborhood_summary(df = annot[, c("fov", "totalcounts", "celltype")],
                                      neighbors = neighbors)
test_that("dataframe_neighborhood_summary returns a list of correct results", {
  expect_true(is.list(tmp))
  expect_true(is.vector(tmp$fov))
  expect_true(is.vector(tmp$totalcounts))
  expect_true(is.matrix(tmp$celltype))
  expect_true(dim(tmp$celltype)[2] == length(unique(annot$celltype)))
})

tmp <- get_neighborhood_expression(counts = counts, neighbors = neighbors)
test_that("get_neighborhood_expression is correct", {
  expect_identical(dim(tmp), dim(counts))
  expect_equal(unname(tmp[2, "SPP1"]),
               sum(counts[neighbors[2, ] > 0, "SPP1"]) / sum(neighbors[2, ] != 0))
})

#### test conditional correlation wrapper and subsidiary functions -----------------------

tmp <- calcSpatialCor (counts = counts,
                  conditionon = annot[, c("fov", "totalcounts", "celltype")],
                  neighbors = NULL,
                  xy = xy,
                  k = 50,
                  radius = NULL,
                  tissue = annot$fov,
                  roundcortozero = 0.1,
                  max_cells = 1e5,
                  verbose = FALSE)
test_that("calcSpatialCor  is correct", {
  expect_true(class(tmp$condcor) == "dsCMatrix")
  expect_true(is.matrix(tmp$env))
  expect_true(class(tmp$neighbors) == "dgCMatrix")
  expect_identical(dim(tmp$condcor), rep(ncol(counts), 2))
  expect_identical(dim(tmp$env), dim(counts))
  expect_identical(dim(tmp$neighbors), rep(nrow(xy), 2))
})

#### test module building wrapper -----------------------------
modules <- defineModules(condcor = tmp$condcor,
                              env = tmp$env,
                              min_module_size = 2,
                              max_module_size = 20,
                              gene_weighting_rule = "inverse_sqrt")
test_that("defineModules is correct", {
  expect_true(is.list(modules))
  expect_true(is.list(modules$modules))
  expect_true(is.character(modules$modules[[1]]))
  expect_true(is.list(modules$weights))
  expect_true(is.numeric(modules$weights[[1]]))
  expect_true(is.data.frame(modules$weightsdf))
})

#### test scoring wrapper --------------------------------

scores <- scoreModules(counts = counts,
                       weights = modules$weights,
                       neighbors = neighbors)
test_that("scoreModules is correct", {
  expect_true(is.list(scores))
  expect_true(is.matrix(scores$scores_env))
  expect_true(is.matrix(scores$scores_sc))
})


#### test cell type attribution scoring ----------------------------

attribution <- cellTypeAttribution(
  modulescores = scores$scores_env,
  weights = modules$weights,
  counts = counts,
  celltype = annot$celltype,
  neighbors = neighbors,
  nsub = 1000,
  verbose = TRUE)

test_that("score_module_celltype_involvement is correct", {
  expect_true(is.list(attribution))
  expect_true(is.matrix(attribution$involvescores))
  expect_identical(dimnames(attribution$involvescores)[[1]], names(modules$weights))
  expect_true(length(setdiff(dimnames(attribution$involvescores)[[2]], unique(annot$celltype))) == 0)
  expect_true(length(setdiff(unique(annot$celltype), dimnames(attribution$involvescores)[[2]])) == 0)
  expect_true(is.matrix(attribution$attributionmats[[1]]))
  expect_identical(colnames(attribution$attributionmats[[1]]), names(modules$weights[[1]]))
  expect_true(length(setdiff(rownames(attribution$attributionmats[[1]]), unique(annot$celltype))) == 0)
  expect_true(length(setdiff(unique(annot$celltype), rownames(attribution$attributionmats[[1]]))) == 0)
})


#### test plots ---------------------------------

plotCorrelationNetwork(x = res$condcor, modules = res$modules, genes = NULL, show_gene_names = FALSE)
plotCorrelationNetwork(x = res$condcor, modules = res$modules, genes = NULL, show_gene_names = TRUE, corthresh = 0.1) 
  