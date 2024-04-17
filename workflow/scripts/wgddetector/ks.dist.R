#!/usr/bin/env Rscript
##  ref. generated 03.hierarchial_clustering/{cluster}/ks.dist.R from WGDdetector, Yang et al. 2019
inputdata <- read.table(snakemake@input[["dist"]])
datadist <- as.dist(inputdata)
datahclust <- hclust(datadist)
out_id <- cutree(datahclust, k = 2)
write.table(t(labels(out_id)), snakemake@output[["clust"]], append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(t(out_id), file = snakemake@output[["clust"]], append = TRUE, row.names = FALSE, col.names = FALSE)