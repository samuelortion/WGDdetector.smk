"""
ref. phase.id.pl

Filter the CDS FASTA file to keep only sequences that appear in the protein FASTA file
"""


rule filter_cds_pep_seqids:
    input:
        cds_fasta=input_cds,
        pep_fasta=input_pep,
    output:
        cds_fasta="{name}.cds.filt.fa",
        pep_fasta="{name}.pep.filt.fa",
        id_table="{name}.filt.id.table",
    conda:
        "../envs/bioperl.yaml"
    log:
        stderr=str(logdir / "filter_cds_pep_seqids" / "{name}.stderr"),
        stdout=str(logdir / "filter_cds_pep_seqids" / "{name}.stdout"),
    shell:
        """
        perl "workflow/scripts/phase.id.pl" "{input.cds_fasta}" "{input.pep_fasta}" "{output.cds_fasta}" "{output.pep_fasta}" "{output.id_table}" > "{log.stdout}" 2> "{log.stderr}"
        """
