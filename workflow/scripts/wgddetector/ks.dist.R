#!/usr/bin/env Rscript
# ref. generated 03.hierarchial_clustering/{cluster}/ks.dist.R from WGDdetector, Yang et al. 2019
# Amend to generate multiple cluster once
inputdata <- read.table(snakemake@input[["dist"]])
datadist <- as.dist(inputdata)
datahclust <- hclust(datadist)

k_start <- 1 # starting cut k value
k_end <- nrow(inputdata) # should be equivalent to `my $linenum=`wc -l $in` in wgddetector/hclust_ks.pl
num <- snakemake@params["num"]
if (num != "all") {
    k_start <- snakemake@wildcards["num"]
    k_start <- as.numeric(k_start)
    if (k_start == NA) {
        message("Error: selected hclust cuttree k is not a number nor \"all\"")
        quit(status=1)
    }
    k_end <- k_start # a single k value is used.
}

for (k in k_start:k_end) {
    out_id <- cutree(datahclust, k=k)
    write.table(t(labels(out_id)), snakemake@output[["clust"]], append = TRUE, row.names = FALSE, col.names = FALSE)
    write.table(t(out_id), file = snakemake@output[["clust"]], append = TRUE, row.names = FALSE, col.names = FALSE)
}