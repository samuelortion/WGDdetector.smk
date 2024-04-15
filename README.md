# WGDdetector.smk
Snakemake port of [WGDdetector](https://github.com/yongzhiyang2012/WGDdetector) pipeline, by Yang _et al._ 2019

WGDDetector is a pipeline for whole genome duplication (WGD) detection with genome or transcriptome annotations.

> [!Note] 
> The goal of this repository is to provide an implementation that matches the results given by the published pipeline while being easier to use; and to learn more doing so.

> [!Warning]
> This repository is a WIP. It currently does not implement all WGDdetector steps.

## Install

Clone this repository (use `--recurse-submodules`)
```bash
git clone --recurse-submodules https://github.com/samuelortion/WGDdetector.smk.git
cd WGDdetector
```

### Dependencies

- Snakemake

If you want to use conda to install Snakemake, we provide a small conda `environment.yaml` file with a minimal set of dependencies.


## Run

Example command line:
```bash
snakemake --cores 1 --snakefile ./workflow/Snakefile --use-conda
```

A small wrapper on the snakemake CLI is provided in `cli.py`, to offer a similar interface and option as the original WGDdetector, with the drawbacks of being less versatile.
Example of cli adapted from `wgddetector/example/00.run.sh`: 
```bash
python3 cli.py --input_cds test.cds.fa --input_pep test.pep.fa --output_dir output --tmp_dir tmp --thread_num 4 --cluster_engine mmseqs2
```

We rather recommend to use the `snakemake` CLI.

## References

> Yang Y, Li Y, Chen Q, Sun Y, Lu Z: WGDdetector: a pipeline for detecting whole genome duplication events using the genome or transcriptome annotations. BMC Bioinformatics 2019, 20(1):75.