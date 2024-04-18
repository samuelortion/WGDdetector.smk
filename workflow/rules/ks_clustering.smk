"""
Ks based clustering of gene families

Steps:
- Split 
- Calculate Ks using Bioperl
- Cluster genes based on Ks
"""

ks_estimate_matrix_folder = (
    outdir
    / config["folder_names"]["ks_estimate_step"]
    / config["folder_names"]["ks_estimate_matrix_step"]
)
ks_cluster_step_folder = outdir / config["folder_names"]["hierachical_clustering_step"]

# TODO
# rule calculate_ks:
#     input:
#         align="{cluster}.cds.fa.align",
#         fasta="{cluster}.cds.fa"
#     output:


rule hclust_ks_dict_all:
    input:
        dist=ks_estimate_matrix_folder / "{cluster}" / "ks.dist",
    output:
        clust=ks_estimate_matrix_folder / "{cluster}" / "ks.dist.hcluster",
    conda:
        "../envs/R.yaml"
    log:
        stdout=logdir / "hclust_ks_dict_all" / "{cluster}.stdout",
        stderr=logdir / "hclust_ks_dict_all" / "{cluster}.stderr",
    params:
        num="all",
    script:
        "../scripts/wgddetector/ks.dist.R"


rule hclust_ks_dict_one:
    input:
        dist=ks_cluster_step_folder / "{cluster}" / "ks.dist",
    output:
        clust=ks_cluster_step_folder / "{cluster}" / "ks.dist.hcluster",
    conda:
        "../envs/R.yaml"
    log:
        stdout=logdir / "hclust_ks_dict_one" / "{cluster}.stdout",
        stderr=logdir / "hclust_ks_dict_one" / "{cluster}.stderr",
    params:
        num="2",  # FIXME: was this hard coded?
    script:
        "../scripts/wgddetector/ks.dist.R"
