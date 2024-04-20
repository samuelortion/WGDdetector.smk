"""
All against all protein local alignment with MMSeqs2
"""


rule mmseqs2_createdb:
    input:
        pep_fasta=mmseqs2_input_fasta,
    params:
        db=mmseqs2_db_path,
    output:
        db_files=multiext(
            str(mmseqs2_db_path),
            "",
            "_h",
            "_h.dbtype",
            "_h.index",
            ".index",
            ".lookup",
            ".source",
        ),
    conda:
        "../envs/mmseqs2.yaml"
    log:
        stderr=logdir / "mmseqs2_createdb.stderr",
        stdout=logdir / "mmseqs2_createdb.stdout",
    shell:
        """
        mmseqs createdb "{input.pep_fasta}" "{params.db}" > "{log.stdout}" 2> "{log.stderr}"
        """


rule mmseqs2_search:
    input:
        pep_fasta=filtered_pep_fasta,
        db_files=multiext(
            str(mmseqs2_db_path),
            "",
            "_h",
            "_h.dbtype",
            "_h.index",
            ".index",
            ".lookup",
            ".source",
        ),
    params:
        db=mmseqs2_db_path,
        out_db=mmseqs2_out_db_path,
        tmp_dir=tmpdir,
    threads: config["threads"]
    output:
        out_db_files=multiext(str(mmseqs2_out_db_path), "", ".dbtype", ".index"),
    conda:
        "../envs/mmseqs2.yaml"
    log:
        stderr=logdir / "mmseqs2_search.stderr",
        stdout=logdir / "mmseqs2_search.stdout",
    shell:
        """
        mmseqs search "{params.db}" "{params.db}" "{params.out_db}" "{params.tmp_dir}" --search-type 2 --threads "{threads}" > "{log.stdout}" 2> "{log.stderr}"
        """


rule mmseqs2_convertalis:
    input:
        db_files=multiext(
            str(mmseqs2_db_path),
            "",
            "_h",
            "_h.dbtype",
            "_h.index",
            ".index",
            ".lookup",
            ".source",
        ),
        out_db_files=multiext(str(mmseqs2_out_db_path), "", ".dbtype", ".index"),
    params:
        db=mmseqs2_db_path,
        out_db=mmseqs2_out_db_path,
    threads: config["threads"]
    output:
        blast_tsv=mmseqs2_blast_tsv,
    conda:
        "../envs/mmseqs2.yaml"
    log:
        stderr=logdir / "mmseqs2_convertalis.stderr",
        stdout=logdir / "mmseqs2_convertalis.stdout",
    shell:
        """
        mmseqs convertalis "{params.db}" "{params.db}" "{params.out_db}" "{output.blast_tsv}" --format-mode 2 --threads "{threads}" > "{log.stdout}" 2> "{log.stderr}"
        """
