"""
Ks based clustering of gene families
"""


rule ks_dist:
    input:
        dist="{cluster}.ks.dist",
    output:
        clust="{cluster}.ks.dist.hcluster",
    conda:
        "../envs/R.yaml"
    log:
        stdout=logdir / "ks_dist" / "{cluster}.stdout",
        stderr=logdir / "ks_dist" / "{cluster}.stderr",
    script:
        "../scripts/wgddetector/ks.dist.R"
