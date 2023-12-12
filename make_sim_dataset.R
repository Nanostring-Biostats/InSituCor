# (execute this within the project directory for this data - SMI-0016)

library(InSituType)

# load data:
clusts = readRDS( "processed_data/clusts.RDS")
cols <- readRDS("processed_data/cols.RDS")
load("processed_data/fundamental data.RData")
rownames(customlocs) = annot$cell_ID
counts <- Matrix::t(raw)



# take just one tissue and 50 genes:
keep = annot$tissue == "SP17_2566"
table(keep)

# keep the top genes per cell type:
genes = c()
profiles = InSituType::Estep(counts = counts, neg = annot$neg, clust = clusts)
for (cell in unique(clusts)) {
  ratios = profiles[, cell] / apply(profiles[, setdiff(colnames(profiles), cell)], 1, max)
  genes = c(genes, names(ratios)[order(ratios, decreasing = TRUE)[1:2]])
}


# Get a subset of data:
simdat <- list()
simdat$counts = counts[keep, genes]
simdat$xy = customlocs[keep, ]
simdat$annot = annot[keep, c("cell_ID", "negmean", "totalcounts", "tissue")]
simdat$celltype = clusts[keep]

# now resample the counts within each cell type:
set.seed(0)
for (cell in unique(simdat$celltype)) {
  inds = which(simdat$celltype == cell)
  simdat$counts[inds, ] = simdat$counts[sample(inds, length(inds), replace = F), ]
}

# now take an FOV and evelate certain genes:
z = simdat$xy[,1] + simdat$xy[,2]
# perturb genes by it:
perturbgenes = genes[1:10]
simdat$perturbgenes = perturbgenes
set.seed(0)
for (gene in perturbgenes) {
  simdat$counts[, gene] = simdat$counts[, gene] + rpois(nrow(simdat$counts), lambda = pmax(z/4, 0))
}
str(simdat)

plot(simdat$xy, pch = 16, cex = 0.5, col = viridis_pal(option="B")(6)[1 + pmin(simdat$counts[, gene], 5)])
save(simdat, file = "simdat.RData")
