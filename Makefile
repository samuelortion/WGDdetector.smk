test:
	snakemake --cores 1 --snakefile ./workflow/Snakefile --rerun-incomplete --use-conda

dag:
	snakemake --dag > tmp/dag.dot
	dot -O tmp/dag.dot -T png
