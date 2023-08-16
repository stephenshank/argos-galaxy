#!/bin/bash

bp_acc=$1
echo Making TSV for $bp_acc...
alias sm='snakemake --rerun-triggers mtime --rerun-incomplete -j 3 -k'
assembly_accessions=data/ncbi/bioproject/$bp_acc/links/assembly/accessions.txt
assembly_content=$(cat $assembly_accessions)
if [ "$assembly_content" != "null" ]; then
  echo "Fetching assembly data for $bp_acc..."
  sm $assembly_accessions
  sm $(sed 's/.*/data\/ncbi\/assembly\/&\/biosample.txt/' $assembly_accessions)
  tsv=data/ncbi/bioproject/$bp_acc/assemblyqc_input.tsv
  sm $tsv
else
  echo "No assemblies for $bp_acc, skipping..."
fi
echo ...made $tsv!
