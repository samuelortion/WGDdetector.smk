"""
ref. phase.id.pl

Filter the CDS FASTA file to keep only sequences that appear in the protein FASTA file
"""


rule filter_cds_pep_seqids:
    input:
        cds_fasta=input_cds,
        pep_fasta=input_pep,
    output:
        cds_fasta=filtered_pep_fasta,
        pep_fasta=filtered_cds_fasta,
        id_table=genenum_seqid_mapping_table,
    conda:
        "../envs/bioperl.yaml"
    log:
        stderr=logdir / "filter_cds_pep_seqids.stderr",
        stdout=logdir / "filter_cds_pep_seqids.stdout",
    shell:
        """
        perl "workflow/scripts/phase.id.pl" "{input.cds_fasta}" "{input.pep_fasta}" "{output.cds_fasta}" "{output.pep_fasta}" "{output.id_table}" > "{log.stdout}" 2> "{log.stderr}"
        """
