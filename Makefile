dag:
	snakemake --dag > tmp/dag.dot
	dot -O tmp/dag.dot -T png
