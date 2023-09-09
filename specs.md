

#### Specs for sparc:
expect_identical(colnames(res$modules), c("module", "gene", "weight"))
  expect_true(nrow(res$scores_env) == nrow(counts))
  expect_true(nrow(res$scores_sc) == nrow(counts))
  expect_true(is(res$condcor, 'sparseMatrix'))
  expect_true(is.matrix(res$celltypeinvolvement))
  expect_true(is.matrix(res$attributionmats[[1]]))
  
  
- Returns a matrix of environment scores aligning to counts -- test: unit_tests.R#L48
- Returns a matrix of single cell scores aligning to counts -- test: unit_tests.R#L49
- Returns a sparse matrix of conditional correlations -- test: unit_tests.R#L50
- Returns a matrix of cell type vs module attribution stats -- test: unit_tests.R#L51
- Returns a list of matrices holding cell type vs gene attribution stats -- test: unit_tests.R#L52




#### Specs for calcSpatialCor:
 - Returns a sparse matrix holding conditional correlations -- test: unit_tests.R#L128,131
 - Returns a matrix of environment expression values -- test: unit_tests.R#L129,132
 - Returns a sparse matrix of neighbor relationships -- test: unit_tests.R#L130,133
 

#### Specs for defineModules:
- Returns a list of module memberships   -- test: unit_tests.R#L144,145
- Returns a list of module weights -- test: unit_tests.R#L146,147
- Returns a data frame summarizing weights and memberships  -- test: unit_tests.R#L148


#### Specs for scoreModules:
- Returns a matrix of environment scores  -- test: unit_tests.R#L196
- Returns a matrix of single cell scores  -- test: unit_tests.R#L197

#### Specs for cellTypeAttribution:
- Returns a matrix of cell type x module attribution scores  -- test: unit_tests.R#L176,177,178,179
- Returns a list of cell type x gene attribution scores  -- test: unit_tests.R#L180,181,182,182


#### Specs for plotCorrelationNetwork
- Draws a plot  -- test: unit_tests.R#L189,190
