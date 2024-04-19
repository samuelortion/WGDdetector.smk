rule blast_blastp_all:
    input:
        query=outdir / "{name}.fa",
        db=expand(
            "{db_path}{ext}",
            db_path=blast_db_path,
            ext=[
                ".pdb",
                ".phr",
                ".pin",
                ".pog",
                ".pos",
                ".pot",
                ".psq",
                ".ptf",
                ".pto",
            ],
        ),
    output:
        tsv=cluster_step_folder / "blastp" / "{name}.blastp.tsv",
    params:
        db=blast_db_path,
        seg="no",
        evalue=config["minimum_homology_evalue"],
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
        makeblastdb -in "{input.fasta}" -parse_seqids -dbtype prot -out "{params.db}" > "{log.stdout}" 2> "{log.stderr}" # makeblastdb returns ERROR even if it ran successfully on some (recent?) versions: exit 0 to avoid stopping the pipeline
        """
