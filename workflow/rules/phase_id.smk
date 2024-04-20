"""
ref. phase.id.pl

Filter the CDS and protein FASTA files to keep only sequences of genes that appear in both files
"""


rule filter_cds_pep_seqids:
    """Phase raw input cds and pep fasta files (keep only matching seq IDs)"""
    input:
        cds_fasta="{name}.cds.fa",
        pep_fasta="{name}.pep.fa",
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
        perl "workflow/scripts/wgddetector/phase.id.pl" "{input.cds_fasta}" "{input.pep_fasta}" "{output.cds_fasta}" "{output.pep_fasta}" "{output.id_table}" > "{log.stdout}" 2> "{log.stderr}"
        """
