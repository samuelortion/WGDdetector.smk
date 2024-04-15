blast_db_path = outdir / "blastp" / "db" / f"{run_name}.pep.filt.fa"


rule blast_blastp_all:
    input:
        multiext(
            str(blast_db_path),
            ".pdb",
            ".phr",
            ".pin",
            ".pjs",
            ".pog",
            ".pos",
            ".pot",
            ".psq",
            ".ptf",
            ".pto",
        ),
        query=outdir / "{name}",
    output:
        tsv=outdir / "blastp/{name}.blastp.tsv",
    params:
        db=blast_db_path,
        seg="no",
        evalue="1e10",  # FIXME: hard-coded in original WGDdetector protein_blastp_cluster.pl
        format="7 std qlen slen",  # we need to keep both lengths
        num_threads=config["threads"],
    conda:
        "../envs/blast.yaml"
    log:
        stderr=logdir / "blast_blastp_all" / "{name}.stdout",
        stdout=logdir / "blast_blastp_all" / "{name}.stderr",
    shell:
        """
        blastp -query "{input.query}" -db "{params.db}" -seg "{params.seg}" -evalue "{params.evalue}" -out "{output.tsv}" -outfmt "{params.format}" -num_threads "{params.num_threads}"  > "{log.stdout}" 2> "{log.stderr}"
        """


rule blast_makeblastdb_protein:
    input:
        fasta=blastp_input_fasta,
    output:
        multiext(
            str(blast_db_path),
            ".pdb",
            ".phr",
            ".pin",
            ".pjs",
            ".pog",
            ".pos",
            ".pot",
            ".psq",
            ".ptf",
            ".pto",
        ),
    params:
        db=blast_db_path,
    conda:
        "../envs/blast.yaml"
    log:
        stderr=logdir / "blast_makeblastdb_protein.stderr",
        stdout=logdir / "blast_makeblastdb_protein.stdout",
    shell:
        """
        makeblastdb -in "{input.fasta}" -parse_seqids -dbtype prot -out "{params.db}" || exit 0  > "{log.stdout}" 2> "{log.stderr}" # makeblastdb returns ERROR even if it ran successfully on some (recent?) versions: exit 0 to avoid stopping the pipeline
        """


rule blast2graphs:
    input:
        blastp_tsv=outdir / "blastp" / "{name}.blastp.tsv",
    output:
        abc_file=outdir / "mcl" / "{name}_nrm_dmls_bit.abc",
    conda:
        "../envs/blast2graphs.yaml"
    log:
        stderr=logdir / "blast2graphs" / "{name}.stderr",
        stdout=logdir / "blast2graphs" / "{name}.stdout",
    shell:
        """
        python3 "workflow/lib/BlastGraphMetrics/blast2graphs.py" "{input.blastp_tsv}" "{wildcards.name}" > "{log.stdout}" 2> "{log.stderr}"
        """
