"""
Ks based clustering of gene families

Steps:
- Split 
- Calculate Ks using Bioperl
- Cluster genes based on Ks

"""

ks_estimate_step_folder = outdir / config["folder_names"]["ks_estimate_step"]
ks_estimate_matrix_folder = (
    outdir / ks_estimate_step_folder / config["folder_names"]["ks_estimate_matrix_step"]
)
ks_cluster_step_folder = outdir / config["folder_names"]["hierarchical_clustering_step"]


checkpoint split_seq:
    """
    Split the input CDS and protein FASTA file according to the previously determined gene clusters
    """
    input:
        cds_file=filtered_cds_fasta,
        pep_file=filtered_pep_fasta,
        cluster_file=cluster_output_file,
    output:
        cluster_cds_fasta=ks_estimate_step_folder / "align" / "{cluster}.input.cds.file",
        cluster_pep_fasta=ks_estimate_step_folder / "align" / "{cluster}.input.pep.file",
    conda:
        "../envs/bioperl.yaml"
    params:
        tmpdir=ks_estimate_step_folder,
    shell:
        """
        perl "workflow/scripts/wgddetector/split_seq.pl" "{input.cds_file}"  "{input.pep_file}" "{input.cluster_file}"  "{params.tmpdir}"
        """


# See cds_muscle, here.


# input function for the rule calculate_ks
def calculate_ks_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.split_seq.get(cluster=wildcards.cluster).output[0].open() as f:
        return ks_estimate_step_folder / "align" / "{cluster}.input.cds.file.align"


rule calculate_ks:
    input:
        align=ks_estimate_step_folder / "align" / "{cluster}.input.cds.file.align",
    output:
        ks=ks_estimate_step_folder / "align" / "{cluster}.input.cds.file.output.ks.gz",
    conda:
        "../envs/bioperl.yaml"
    log:
        stdout=logdir / "calculate_ks" / "{cluster}.stdout",
        stderr=logdir / "calculate_ks" / "{cluster}.stderr",
    params:
        default_ks=10,  # TODO: really use this in the perl script
    shell:
        """
        perl "workflow/scripts/wgddetector/calculate_ks.pl" "{input.align}" "{output.ks}"
        """


# print(ks_estimate_matrix_folder / "{cluster}" / "ks.dist")
# print(ks_estimate_step_folder
#         / "align"
#         / "{cluster}.input.cds.file.output.ks.gz",)
# print(ks_estimate_matrix_folder / "align" / "{cluster}.input.cds.file.output.ks.gz",)


rule ks_dist:
    input:
        ks_list=ks_estimate_step_folder
        / "align"
        / "{cluster}.input.cds.file.output.ks.gz",
    output:
        ks_dist=ks_estimate_matrix_folder / "{cluster}" / "ks.dist",
    shell:
        """
        perl "workflow/scripts/wgddetector/collect_ks.pl" "{input.ks_list}" "{output.ks_dist}"
        """


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
