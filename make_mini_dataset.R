# (execute this within the project directory for this data - SMI-0016)
# load data:
clusts = readRDS( "processed_data/clusts.RDS")
cols <- readRDS("processed_data/cols.RDS")
load("processed_data/fundamental data.RData")
rownames(customlocs) = annot$cell_ID
counts <- Matrix::t(raw)

# what to subsample:
xlim = c(2.1, 3.05)
ylim = c(-8.25, -6.95)
genes = c("ITGAV", "ITGA3", "SPOCK2", "SPP1", "SAT1", "SERPINA1",
          "CALB1",  "IL18", "KIT", "KRT17", "B2M")

keep = ((viz[, 1] < xlim[2]) & (viz[, 1] > xlim[1])) & ((viz[, 2] < ylim[2]) & (viz[, 2] > ylim[1]))

out = list(counts = counts[keep, genes],
           annot = annot[keep, c("cell_ID", "fov", "totalcounts", "negmean")],
           celltype = clusts[keep],
           xy = viz[keep, ])
saveRDS(out, file = "CosMx kidney 2 FOVs 12 genes.RDS")
