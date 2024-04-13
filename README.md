# WGDdetector.smk
Snakemake port of [WGDdetector](https://github.com/yongzhiyang2012/WGDdetector) pipeline, by Yang _et al._ 2019

WGDDetector is a pipeline for whole genome duplication (WGD) detection with genome or transcriptome annotations.

> [!Note] 
> The goal of this repository is to provide an implementation that matches the results given by the published pipeline while being easier to use; and to learn more doing so.

## Install

### Dependencies

- Snakemake

If you want to use conda to install Snakemake, we provide a small conda `environment.yaml` file with a minimal set of dependencies.


## Run

Example command line:
```bash
snakemake --cores 1 --snakefile ./workflow/Snakefile --use-conda
```

## References

> Yang Y, Li Y, Chen Q, Sun Y, Lu Z: WGDdetector: a pipeline for detecting whole genome duplication events using the genome or transcriptome annotations. BMC Bioinformatics 2019, 20(1):75.