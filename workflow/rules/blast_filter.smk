"""
Filter blast output

(Rewriting of a section the protein_mmseqs2_cluster.pl script, in AWK)
"""


rule filter_blast_hsp_evalue:
    """Filter out blast HSPs above a given value"""
    input:
        blast_tsv_raw=blast_tsv_raw,
    output:
        blast_tsv_filtered=blast_tsv_filtered,
    params:
        threshold="1e-10",  # FIXME: handle this hard-coded value properly
    # conda: # TODO: Should I add a awk dependencies yaml file?
    log:
        stderr=logdir / "filter_blast_hsp_evalue.stderr",
        # stdout=logdir / "filter_blast_hsp_evalue.stdout", # stdout is already used to redirect awk output
    shell:
        """
        awk -f "workflow/scripts/filter_blast_hsp_evalue.awk" -v threshold="{params.threshold}" "{input.blast_tsv_raw}" > "{output.blast_tsv_filtered}" 2> "{log.stderr}"
        """
