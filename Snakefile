import xmltodict


wildcard_constraints:
  db="[^/]+",
  id_="[^/]+",
  linked_db="[^/]+"

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

rule assembly_refseq_links:
  input:
    rules.entrez_elink_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/links/{linked_db}/refseq.txt'
  shell:
    '''
      jq -r '.DocumentSummarySet.DocumentSummary | map(select(.SourceDb == "refseq")) | .[].Caption' {input} > {output}
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

rule bioproject_assembly_biosample_links:
  input:
    rules.entrez_elink_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/links/{linked_db}/biosample.txt'
  shell:
    'jq -r ".DocumentSummarySet | .[].BioSampleAccn" {input} > {output}'

rule biosample_sra_links:
  input:
    rules.entrez_elink_summary.output.json
  output:
    'data/ncbi/{db}/{id_}/links/{linked_db}/sra_links.txt'
  shell:
    'sed "s/\@//g" {input} | jq -r ".DocumentSummarySet[].Runs.Run.acc" > {output}'

rule all:
  input:
    # Original ARGOS BioProject assemblies
    'data/ncbi/bioproject/PRJNA231221/links/assembly/summary.json'
