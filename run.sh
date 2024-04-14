#!/usr/bin/env bash
snakemake --cores 1 --snakefile ./workflow/Snakefile --use-conda --rerun-incomplete
