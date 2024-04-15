"""
Filter blast output

(Rewriting of a section the protein_mmseqs2_cluster.pl script, in AWK)
"""


rule filter_evalue_blast_hsp:
    """Filter out blast HSPs above a given value"""
    input:
        "{name}.blast",
    output:
        "{name}_fevalue.blast",
    params:
        threshold="1e-10",  # FIXME: handle this hard-coded value properly
    # conda: # TODO: Should I add a awk dependencies yaml file?
    log:
        stderr=logdir / "filter_blast_hsp_evalue.stderr",
        # stdout=logdir / "filter_blast_hsp_evalue.stdout", # stdout is already used to redirect awk output
    shell:
        """
        awk -f "workflow/scripts/filter_blast_hsp_evalue.awk" -v threshold="{params.threshold}" "{input}" > "{output}" 2> "{log.stderr}"
        """
