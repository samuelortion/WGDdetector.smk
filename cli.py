#!/usr/bin/env python3
"""
A small wrapper around the Snakefile to match wgddetector.pl CLI

ref. https://charlesreid1.github.io/building-snakemake-command-line-wrappers-for-workflows.html
"""

import os
import argparse
import subprocess, sys


thisdir = os.path.abspath(os.path.dirname(__file__))

def config_to_args(config):
    """Convert config dict into a string of --config KEY=VAL"""
    return " ".join(
        [""] + [f"--config {key}={value}" for key, value in config.items()]
    )

def args_to_config(args):
    config = {}
    config["input_cds"] = args.input_cds
    config["input_pep"] = args.input_pep
    config["outdir"] = args.output_dir
    config["tmpdir"] = args.tmp_dir
    config["threads"] = args.thread_num
    config["cluster_engine"] = args.cluster_engine
    return config

def snakemake_cli(snakefile, cores, dryrun, configfile, configoptions):
    """Build the snakemake system call"""
    dryoption = "--dry-run " if  dryrun else ""
    return f"""snakemake --cores "{cores}" --snakefile "{snakefile}" {dryoption} --use-conda \
        --configfile "{configfile}" {configoptions}
    """

def launch(args):
    config = args_to_config(args)
    configoptions = config_to_args(config)

    snakefile = os.path.join(thisdir, "workflow", "Snakefile")
    command = snakemake_cli(
        snakefile,
        args.cores,
        configfile=args.configfile,
        dryrun=args.dry_run,
        configoptions=configoptions,
    )

    return subprocess.run(command, shell = True)


def main():
    parser = argparse.ArgumentParser(
        prog="WGDdetector.smk",
        description="Whole Genome Duplication detection pipeline ported to Snakemake",
        epilog="Sequence IDs within CDS and protein files must be the same.",
    )
    # wgddetector.pl arguments
    parser.add_argument("--input_cds", default=None, help="input CDS in FASTA format")
    parser.add_argument(
        "--input_pep", default=None, help="input protein in FASTA format"
    )
    parser.add_argument("--output_dir", default=None, help="output directory")
    parser.add_argument(
        "--tmp_dir", default=None, help="temporary files directory"
    )  # TODO: decide whether to use this.
    parser.add_argument("--thread_num", default=None, help="number of threads to use")
    parser.add_argument(
        "--cluster_engine",
        choices=["blastp", "mmseqs2"],
        default=None,
        help="local alignment method",
    )

    # Snakemake arguments
    parser.add_argument(
        "--configfile", default="config/config.yaml", help="config file"
    )
    parser.add_argument("-n", "--dry-run")
    parser.add_argument("--cores", default=1)
    
    args = parser.parse_args()

    launch(args)


if __name__ == "__main__":
    main()
