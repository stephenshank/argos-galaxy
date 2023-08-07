import json

import xmltodict


wildcard_constraints:
  db="[^/]+",
  id_="[^/]+",
  linked_db="[^/]+"


def read_lines(input_filepath):
  with open(input_filepath) as input_file:
    lines = [line.strip() for line in input_file.readlines()]
  return lines


rule entrez_fetch_xml:
  output:
    xml='data/ncbi/{db}/{id_}/fetch.xml',
    json='data/ncbi/{db}/{id_}/fetch.json'
  shell:
    '''
      efetch -db {wildcards.db} -id {wildcards.id_} -format xml > {output.xml}
      xmltojson {output.xml} {output.json}
    '''

rule entrez_fetch_fasta:
  output:
    'data/ncbi/nuccore/{id_}/fetch.fasta'
  shell:
    'efetch -db nuccore -id {wildcards.id_} -format fasta> {output}'

rule fasta_header:
  input:
    rules.entrez_fetch_fasta.output[0]
  output:
    'data/ncbi/nuccore/{id_}/header.txt'
  shell:
    'head -n 1 {input} | cut -c 2- > {output}'

rule entrez_summary:
  output:
    xml='data/ncbi/{db}/{id_}/summary.xml',
    json='data/ncbi/{db}/{id_}/summary.json'
  shell:
    '''
      esummary -db {wildcards.db} -id {wildcards.id_} > {output.xml}
      xmltojson {output.xml} {output.json}
    '''

rule summary_biosample_accession:
  input:
    rules.entrez_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/biosample.txt'
  shell:
    'jq -r ".DocumentSummarySet.DocumentSummary.BioSampleAccn" {input} > {output}'

rule assemblyqc_assembly_jq:
  input:
    rules.entrez_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/assemblyqc.json'
  shell:
    '''
      jq '.DocumentSummarySet.DocumentSummary |
        {{
          "organism_name": .SpeciesName,
          "taxonomy_id": .Taxid,
          "assembled_genome_acc": .AssemblyAccession,
          "assembly_file_source": .FtpPath_GenBank
        }}' {input} > {output}
    '''

checkpoint entrez_elink_summary:
  output:
    xml='data/ncbi/{db}/{id_}/links/{linked_db}/summary.xml',
    json='data/ncbi/{db}/{id_}/links/{linked_db}/summary.json'
  params:
    temp_xml='data/ncbi/{db}/{id_}/links/{linked_db}/temp_summary.xml',
    temp_json='data/ncbi/{db}/{id_}/links/{linked_db}/temp_summary.json'
  shell:
    '''
      elink -db {wildcards.db} -id {wildcards.id_} -target {wildcards.linked_db} | esummary > {output.xml}
      dss_xml_cleaner -i {output.xml} -o {params.temp_xml}
      xml_to_json -i {params.temp_xml} -o {params.temp_json}
      dss_json_cleaner -i {params.temp_json} -o {output.json}
    '''


def refseq_inputs(wildcards):
    base = 'data/ncbi/%s/%s/links/nuccore/summary.json'
    variables = (wildcards.db, wildcards.id_)
    return base % variables


rule refseq_linked_accessions:
  input:
    refseq_inputs
  output:
    'data/ncbi/{db}/{id_}/links/nuccore/refseq_accessions.txt'
  shell:
    '''
      jq -r '.DocumentSummarySet | map(select(.SourceDb == "refseq")) | .[].Caption' {input} > {output}
    '''


def sra_inputs(wildcards):
    base = 'data/ncbi/%s/%s/links/sra/summary.json'
    variables = (wildcards.db, wildcards.id_)
    return base % variables


rule sra_linked_run_accessions:
  input:
    sra_inputs
  output:
    'data/ncbi/{db}/{id_}/links/sra/sra_run_accessions.txt'
  shell:
    '''
      jq -r '.DocumentSummarySet[].Runs.Run."@acc"' {input} > {output}
    '''

rule bioproject_assembly_accessions:
  input:
    rules.entrez_elink_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/links/{linked_db}/accessions.txt'
  shell:
    '''
      jq -r '.DocumentSummarySet[].AssemblyAccession' {input} > {output}
    '''

rule biosample_sra_links:
  input:
    rules.entrez_elink_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/links/{linked_db}/sra_links.txt'
  shell:
    'sed "s/\@//g" {input} | jq -r ".DocumentSummarySet[].Runs.Run.acc" > {output}'

rule argos_all_biosample_data:
  input:
    refseq=rules.refseq_linked_accessions.output[0],
    sra=rules.sra_linked_run_accessions.output[0]
  output:
    "data/ncbi/{db}/{id_}/argos_biosample.json"
  run:
    refseq_accessions = read_lines(input.refseq)
    sra_accessions = read_lines(input.sra)
    with open(output[0], 'w') as json_file:
      json.dump({
        'refseq_accessions': refseq_accessions,
        'sra_accessions': sra_accessions
      }, json_file, indent=2)

def aaad_input(wildcards):
    bs_path = "data/ncbi/assembly/%s/biosample.txt" % wildcards.id_
    with open(bs_path) as f:
      biosample_accession = f.read().strip()
    return "data/ncbi/biosample/%s/argos_biosample.json" % biosample_accession

rule argos_all_assembly_data:
  input:
    aaad_input
  output:
    "data/ncbi/{db}/{id_}/argos_assembly.json"
  run:
    with open(input[0]) as json_file:
      biosample_data = json.load(json_file)
    with open(output[0], 'w') as json_file:
      json.dump({
        'assembly_accession': wildcards.id_,
        **biosample_data,
      }, json_file, indent=2)

rule all:
  input:
    # Original ARGOS BioProject assemblies
    'data/ncbi/bioproject/PRJNA231221/links/assembly/summary.json'
