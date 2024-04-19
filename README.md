# WGDdetector.smk
Snakemake port of [WGDdetector](https://github.com/yongzhiyang2012/WGDdetector) pipeline, by Yang _et al._ 2019

WGDDetector is a pipeline for whole genome duplication (WGD) detection with genome or transcriptome annotations.

> [!Note] 
> The goal of this repository is to provide an implementation that matches the results given by the published pipeline while being easier to use; and to learn more doing so.

> [!Warning]
> This repository is a WIP. It currently does not implement all WGDdetector steps.

## Workflow

### Input

- proteome (CDS and protein sequences in FASTA format)

### Output

- gene families (orthomcl format)
- plot of the distribution of pairwise paralogs synonymous substitutions per synonymous site (dS)


### Steps

1. gene similarity search on proteome (BLASTp or MMSeqS2)
2. gene families identification: protein-protein homology graph clustering (MCL)
3. pairwise alignement on CDS (MUSCLE and PAL2NAL)
4. computation of dS
5. identification of sub gene families based on dS (hierarchical clustering with R)
6. sub gene families phasing


## Install

Clone this repository (use `--recurse-submodules`)
```bash
git clone --recurse-submodules https://github.com/samuelortion/WGDdetector.smk.git
cd WGDdetector
```

### Dependencies

- Snakemake

If you want to use conda to install Snakemake, we provide a small conda `environment.yaml` file with a minimal set of dependencies.

If you do not want to use conda, you will have to install the following dependencies on your own:
- perl >5.0, with modules Bioperl and threads
- python
    - [blast2graphs.py](https://github.com/trgibbons/BlastGraphMetrics/blob/master/blast2graphs.py) (It is included as a git submodule, and requires networkx=2.0)
- R >=3.5
- mmseqs2
- ncbi-blast >=2.2.28
- mcl
- muscle
- fastme

<!--TODO: -->
> [!Note]
> I did not check If it works well with all these installed.

## Run

Example command line:
```bash
snakemake --cores 1 --snakefile ./workflow/Snakefile --use-conda
```

A small wrapper on the snakemake CLI is provided in `cli.py`, to offer a similar interface and option as the original WGDdetector. However, we recommend using the `snakemake` CLI instead.

Example of cli adapted from `wgddetector/example/00.run.sh`: 
```bash
python3 cli.py --input_cds test.cds.fa --input_pep test.pep.fa --output_dir output --tmp_dir tmp --thread_num 4 --cluster_engine mmseqs2
```



## References

> Yang Y, Li Y, Chen Q, Sun Y, Lu Z: WGDdetector: a pipeline for detecting whole genome duplication events using the genome or transcriptome annotations. BMC Bioinformatics 2019, 20(1):75.