"""
MUSCLE multiple sequence alignment
"""


rule muscle_align:
    input:
        fasta="{name}.pep.file",
    output:
        align="{name}.pep.file.align",
    conda:
        "../envs/muscle.yaml"
    shell:
        """
        muscle -align "{input.fasta}" -output "{output.align}"
        """


rule pal2nal:
    """Convert protein sequence alignment into the corresponding codon alignment"""
    input:
        palign="{name}.pep.file.align",
        cds="{name}.cds.file",
    output:
        nalign="{name}.cds.file.align",
    conda:
        "../envs/pal2nal.yaml"
    shell:
        """
        pal2nal.pl "{input.palign}" "{input.cds}" -output fasta > "{output.nalign}"
        """
