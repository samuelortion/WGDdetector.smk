"""
Run Markov Clustering (MCL) on an ABC file representing an homology links graph between genes
"""


rule blast2graphs:
    input:
        blastp_tsv=cluster_step_folder / protein_similarity_engine / "{name}.blastp.tsv",
    output:
        abc=cluster_step_folder / "mcl" / "{name}.abc",
    params:
        prefix=lambda wildcards: cluster_step_folder / "mcl" / wildcards.name,
    conda:
        "../envs/blast2graphs.yaml"
    log:
        stderr=logdir / "blast2graphs" / "{name}.stderr",
        stdout=logdir / "blast2graphs" / "{name}.stdout",
    shell:
        """
        python3 "workflow/lib/BlastGraphMetrics/blast2graphs.py" "{input.blastp_tsv}" "{params.prefix}" > "{log.stdout}" 2> "{log.stderr}"
        mv "{params.prefix}_nrm_dmls_bit.abc" "{output.abc}"
        """


rule mcl_genes_abc:
    input:
        abc="{name}.abc",
    output:
        mcl="{name}.mcl",
    conda:
        "../envs/mcl.yaml"
    log:
        stdout=logdir / "mcl_genes_abc" / "{name}.stdout",
        stderr=logdir / "mcl_genes_abc" / "{name}.stderr",  # TODO We could add a thread parameter to speed up mcl.
    shell:
        """
        mcl "{input.abc}" --abc -I 1.5 -o "{output.mcl}" > "{log.stdout}" 2> "{log.stderr}"
        """


rule mcl_to_orthomcl:
    input:
        mcl="{name}.mcl",
    output:
        orthomcl="{name}.orthomcl",
    conda:
        "../envs/bioperl.yaml"
    log:
        stderr=logdir / "mcl_output_to_orthomcl" / "{name}.stderr",
        stdout=logdir / "mcl_output_to_orthomcl" / "{name}.stdout",
    shell:
        """
        perl "workflow/scripts/wgddetector/phase.mcloutp2orthomcl.format.pl" "{input.mcl}" "{output.orthomcl}" > "{log.stdout}" 2> "{log.stderr}"
        """
