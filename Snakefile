import xmltodict


wildcard_constraints:
  db="[^/]+",
  id_="[^/]+",
  linked_db="[^/]+"

rule entrez_summary:
  output:
    xml='data/ncbi/{db}/{id_}/summary.xml',
    json='data/ncbi/{db}/{id_}/summary.json'
  shell:
    '''
      esummary -db {wildcards.db} -id {wildcards.id_} > {output.xml}
      xmltojson {output.xml} {output.json}
    '''

rule entrez_elink_summary:
  output:
    xml='data/ncbi/{db}/{id_}/links/{linked_db}/summary.xml',
    json='data/ncbi/{db}/{id_}/links/{linked_db}/summary.json'
  shell:
    '''
      elink -db {wildcards.db} -id {wildcards.id_} -target {wildcards.linked_db} | esummary > {output.xml}
      xmltojson {output.xml} {output.json}
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

rule all:
  input:
    # Original ARGOS BioProject assemblies
    'data/ncbi/bioproject/PRJNA231221/links/assembly/summary.xml'
