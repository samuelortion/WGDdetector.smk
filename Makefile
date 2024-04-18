test:
	snakemake --cores 1 --snakefile ./workflow/Snakefile --rerun-incomplete --use-conda

dag:
	snakemake --config outdir=./none --dag > tmp/dag.dot
	dot -O tmp/dag.dot -T png
