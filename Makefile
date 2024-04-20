test:
	snakemake --cores 1 --snakefile ./workflow/Snakefile --rerun-incomplete --use-conda

dag:
	snakemake --dag --config outdir=./none > tmp/dag.dot
	dot -O tmp/dag.dot -T png

time:
	time snakemake --cores 1 --snakefile ./workflow/Snakefile --use-conda --config outdir=tmp/results