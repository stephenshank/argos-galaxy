#!/bin/bash

bp_acc=$1
echo Making TSV for $bp_acc...
alias sm='snakemake --rerun-triggers mtime --rerun-incomplete -j 1 -k'
assembly_accessions=data/ncbi/bioproject/$bp_acc/links/assembly/accessions.txt
sm $assembly_accessions
sm $(sed 's/.*/data\/ncbi\/assembly\/&\/biosample.txt/' $assembly_accessions)
tsv=data/ncbi/bioproject/$bp_acc/assemblyqc_input.tsv
sm $tsv
echo ...made $tsv!
