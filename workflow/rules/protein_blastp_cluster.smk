rule blast_blastp_all:
    input:
        expand(
            "{db}{ext}",
            db=blast_db_path,
            ext=[
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
            ],
        ),
        query=proteome_fasta,
    output:
        csv=blastp_csv,
    params:
        db=blast_db_path,
        seg="yes",
        evalue="1e10",  # FIXME: hard-coded in original WGDdetector protein_blastp_cluster.pl
        format="7 std qlen slen",  # we need to keep both lengths
        num_threads=config["threads"],
    conda:
        "../envs/blast.yaml"
    log:
        stderr=logdir / "blast_blastp_all.stderr",
        stdout=logdir / "blast_blastp_all.stdout",
    shell:
        """
        blastp -query "{input.query}" -db "{params.db}" -seg "{params.seg}" -evalue "{params.evalue}" -out "{output.csv}" -outfmt "{params.format}" -num_threads "{params.num_threads}"  > "{log.stdout}" 2> "{log.stderr}"
        """


rule blast_makeblastdb_protein:
    input:
        fasta=proteome_fasta,
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
