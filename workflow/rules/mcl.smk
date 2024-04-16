"""
Run Markov Clustering (MCL) on an ABC file representing an homology links graph between genes
"""


rule mcl_genes_abc:
    input:
        "{name}.abc",
    output:
        "{name}.mcl",
    conda:
        "../envs/mcl.yaml"
    log:
        stdout=logdir / "mcl_genes_abc" / "{name}.stdout",
        stderr=logdir / "mcl_genes_abc" / "{name}.stderr",
    shell:
        """
        mcl "{input}" --abc -I 1.5 -o "{output}" > "{log.stdout}" 2> "{log.stderr}"
        """


rule mcl_output_to_orthomcl:
    input:
        "{name}.mcl",
    output:
        "{name}.orthomcl",
    conda:
        "../envs/bioperl.yaml"
    log:
        stderr=logdir / "mcl_output_to_orthomcl" / "{name}.stderr",
        stdout=logdir / "mcl_output_to_orthomcl" / "{name}.stdout",
    shell:
        """
        perl "workflow/scripts/phase.mcloutp2orthomcl.format.pl" "{input}" "{output}" > "{log.stdout}" 2> "{log.stderr}"
        """
