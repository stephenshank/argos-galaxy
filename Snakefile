import json
import csv

import xmltodict


wildcard_constraints:
  db="[^/]+",
  id_="[^/]+",
  linked_db="[^/]+"


def read_lines(input_filepath):
  with open(input_filepath) as input_file:
    lines = [line.strip() for line in input_file.readlines()]
  return lines


def read_json(input_filepath):
  with open(input_filepath) as input_file:
    d = json.load(input_file)
  return d


def write_json(the_json, filepath):
  with open(filepath, 'w') as json_file:
    json.dump(the_json, json_file)


rule taxonomy_bioprojects:
  output:
    xml="data/ncbi/{db}/{id_}/tax_bps.xml",
    json="data/ncbi/{db}/{id_}/tax_bps.json",
    accessions="data/ncbi/{db}/{id_}/tax_bps.txt"
  params:
    temp_xml='data/ncbi/{db}/{id_}/temp_tax_bps.xml',
    temp_json='data/ncbi/{db}/{id_}/temp_tax_bps.json'
  shell:
    '''
      esearch -db bioproject -query "txid{wildcards.id_}[ORGN]" | esummary > {output.xml}
      dss_xml_cleaner -i {output.xml} -o {params.temp_xml}
      xml_to_json -i {params.temp_xml} -o {params.temp_json}
      dss_json_cleaner -i {params.temp_json} -o {output.json}
      jq -r '.DocumentSummarySet[].Project_Acc' {output.json} > {output.accessions}
    '''

rule taxonomy_bioprojects_with_assemblies:
  input:
    rules.taxonomy_bioprojects.output.accessions
  output:
    "data/ncbi/{db}/{id_}/assembly_tax_bps.txt"
  run:
    with open(input[0]) as tax_bp_file:
      tax_bps = [l.strip() for l in tax_bp_file.readlines()]
    bps_with_assemblies = []
    for tax_bp in tax_bps:
      base = 'data/ncbi/bioproject/%s/links/assembly/accessions.txt'
      with open(base % tax_bp) as links_file:
        assemblies = links_file.read()
      if assemblies != 'null\n':
        bps_with_assemblies.append(tax_bp)
    with open(output[0], 'w') as output_file:
      output_file.write('\n'.join(bps_with_assemblies))

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

rule entrez_elink_summary:
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
      jq -r '
        .DocumentSummarySet |
        if type == "array" then
          map(select(.SourceDb == "refseq" // .SourceDb == "insd")) |
          .[].Caption
        else
          .Caption
        end' {input} > {output}
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
      jq -r '.DocumentSummarySet |
        if type == "array" then
          .[].Runs.Run |
            if type == "array" then
              .[]."@acc"
            else
              ."@acc"
            end
        else
          .Runs.Run |
            if type == "array" then
              .[]."@acc"
            else
              ."@acc"
            end
        end' {input} > {output}
    '''

rule count_sra_run_accessions_json:
  input:
    rules.sra_linked_run_accessions.output[0]
  output:
    'data/ncbi/{db}/{id_}/sra_run_counts.json'
  run:
    lines = [
      line for line in read_lines(input[0]) if line != 'null'
    ]
    with open(output[0], 'w') as output_file:
      json.dump({
        wildcards.id_: len(lines)
      }, output_file)

def count_sra_runs_in_assembly_bps_for_taxon_input(wildcards):
  base = 'data/ncbi/taxonomy/%s/assembly_tax_bps.txt'
  bp_accessions = read_lines(base % wildcards.id_)
  return expand(
    'data/ncbi/bioproject/{bp_accession}/sra_run_counts.json',
    bp_accession=bp_accessions
  )

rule count_sra_runs_in_assembly_bps_for_taxon:
  input:
    count_sra_runs_in_assembly_bps_for_taxon_input
  output:
    'data/ncbi/{db}/{id_}/sra_run_counts_for_assembly_bps.json'
  run:
    all_sra_run_counts = {}
    for input_filepath in input:
      sra_run_counts = read_json(input_filepath)
      all_sra_run_counts.update(sra_run_counts)
    write_json(all_sra_run_counts, output[0])

rule bioproject_assembly_accessions:
  input:
    rules.entrez_elink_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/links/{linked_db}/accessions.txt'
  shell:
    '''
      jq -r ' .DocumentSummarySet |
        if type == "array" then
          .[].AssemblyAccession
        else
          .AssemblyAccession
        end
      ' {input} > {output}
    '''

rule bioproject_biosample_accessions:
  input:
    rules.entrez_elink_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/links/{linked_db}/biosample_accessions.txt'
  shell:
    'jq -r ".DocumentSummarySet[].Accession" {input} > {output}'

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

def asai_input(wildcards):
    with open('data/ncbi/bioproject/%s/links/assembly/accessions.txt' % wildcards.id_) as ids_file:
      accessions = [l.strip() for l in ids_file.readlines()]
    return expand(
      'data/ncbi/assembly/{id_}/argos_assembly.json',
      id_=accessions
    )

rule argos_style_assemblyqc_input:
  input:
    asai_input
  output:
    'data/ncbi/{db}/{id_}/assemblyqc_input.tsv'
  run:
    csv_file = open(output[0], 'w')
    csv_writer = csv.DictWriter(
      csv_file,
      fieldnames=[
        'assembly_acc', 'sra_acc', 'refseq_acc'
      ],
      delimiter="\t"
    )
    csv_writer.writeheader()
    for assembly_metadata_filepath in input:
      with open(assembly_metadata_filepath) as json_file:
        assembly_metadata = json.load(json_file)
      csv_writer.writerow({
        'assembly_acc': assembly_metadata['assembly_accession'],
        'refseq_acc': ','.join(assembly_metadata['refseq_accessions']),
        'sra_acc': ','.join(assembly_metadata['sra_accessions'])
      })
    csv_file.close()

rule all:
  input:
    # Original ARGOS BioProject assemblies
    'data/ncbi/bioproject/PRJNA231221/links/assembly/summary.json'
