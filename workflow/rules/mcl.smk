"""
Run Markov Clustering (MCL) on an ABC file representing an homology links graph between genes
"""


rule mcl_genes_abc:
    input:
        abc_file=abc_file,
    output:
        mcl_file=mcl_file,
    conda:
        "../envs/mcl.yaml"
    log:
        stderr=logdir / "mcl_genes_abcs.stderr",
        stdout=logdir / "mcl_genes_abc.stdout",
    shell:
        """
        mcl "{input.abc_file}" --abc -I 1.5 -o "{output.mcl_file}" > "{log.stdout}" 2> "{log.stderr}"
        """


rule mcl_output_to_orthomcl:
    input:
        mcl_file=mcl_file,
    output:
        orthomcl_file=orthomcl_file,
    conda:
        "../envs/bioperl.yaml"
    log:
        stderr=logdir / "mcl_output_to_orthomcl.stderr",
        stdout=logdir / "mcl_output_to_orthomcl.stdout",
    shell:
        """
        perl "workflow/scripts/phase.mcloutp2orthomcl.format.pl" "{input.mcl_file}" "{output.orthomcl_file}" > "{log.stdout}" 2> "{log.stderr}"
        """
