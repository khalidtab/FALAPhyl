scattergather:
    bootstraps=config["sparcc_bootstrap"][0]

checkpoint makeCore: # SparCC works better when the feature table has less zeros. Therefore, this script will filter the feature tables to the desired level of "coreness" based on the input file values
   version: "1.0"
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      biom="data/biom/{sample}_temp.biom",
      map="data/{sample}.txt"
   params:
      threshold=config["threshold"][0]
   output:
      directory("data/network/{sample}–{group}/core/")
   message: "Making core files for {wildcards.sample} {wildcards.group}"
   shell:
      "mkdir -p {output} && "
      "Rscript --vanilla ./workflow/scripts/make_core.R -i {input.biom} -m {input.map} -c {wildcards.group} -t {params.threshold} -o {output}/ "

rule fastspar: # SparCC (through the fastpar algorithm) will be done on the core feature table. This is an intermediate step that will solve the output file name based on the definition below, and as such solves the input based on the output
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      "data/network/{sample}–{group}/core/{groupcategory}.tsv"
   output:
      "data/network/{sample}–{group}/corr/{groupcategory}.tsv"
   params:
      sparcc=config["sparcc_corr"][0]
   message: "Calculating SparCC on core of {wildcards.sample} {wildcards.group} {wildcards.groupcategory}"
   shell:
      "fastspar -c {input} -r {output} -a {output}.cov.tsv -e {params.sparcc} -y >/dev/null && "
      "rm {output}.cov.tsv"

rule bootstrap: # To calculate P-values, bootstraps need to be done from the core feature table. This script accomplishes this.
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      "data/network/{sample}–{group}/core/{groupcategory}.tsv"
   output:
      temporary("data/network/{sample}–{group}/bootstrap_corr_{groupcategory}_done.txt")
   params:
      bootstrap=config["sparcc_bootstrap"][0]
   message: "Making bootstrap on {wildcards.sample} {wildcards.group} {wildcards.groupcategory}"
   shell:
      '''
      mkdir -p data/network/{wildcards.sample}–{wildcards.group}/bootstrap–{wildcards.groupcategory}/fake_data &&
      fastspar_bootstrap --otu_table {input} --number {params.bootstrap} --prefix data/network/{wildcards.sample}–{wildcards.group}/bootstrap–{wildcards.groupcategory}/fake_data/ >/dev/null 2>/dev/null &&
      ls data/network/{wildcards.sample}–{wildcards.group}/bootstrap–{wildcards.groupcategory}/fake_data/* | parallel fastspar --otu_table {{}} --correlation data/network/{wildcards.sample}–{wildcards.group}/bootstrap–{wildcards.groupcategory}/fake_data/cor_{{/}} --covariance data/network/{wildcards.sample}–{wildcards.group}/bootstrap–{wildcards.groupcategory}/fake_data/cov_{{/}} -i 5 -y >/dev/null 2>/dev/null &&
      rm data/network/{wildcards.sample}–{wildcards.group}/bootstrap–{wildcards.groupcategory}/fake_data/cov_* data/network/{wildcards.sample}–{wildcards.group}/bootstrap–{wildcards.groupcategory}/fake_data/_* &&
      touch {output}
      '''



rule pvalues: # Calculates p-values by using the correlations from the core feature table, and the boostraps
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      bootstrapdone="data/network/{sample}–{group}/bootstrap_corr_{groupcategory}_done.txt",
      core="data/network/{sample}–{group}/core/{groupcategory}.tsv",
      corr="data/network/{sample}–{group}/corr/{groupcategory}.tsv"
   output:
      pvalue="data/network/{sample}–{group}–corr/pvalue–{groupcategory}.tsv"
   log: 
      "data/logs/network_pvalue_{sample}–{group}–{groupcategory}.txt"
   params:
      bootstrap=config["sparcc_bootstrap"][0]
   message: "pvalues calculation for the features in networks of {wildcards.sample} {wildcards.group} {wildcards.groupcategory}"
   shell:
      '''
      fastspar_pvalues --otu_table {input.core} --correlation {input.corr} --prefix data/network/{wildcards.sample}–{wildcards.group}/bootstrap–{wildcards.groupcategory}/fake_data/cor__ --permutations {params.bootstrap} --outfile {output.pvalue} > {log}
      '''


rule nodes_and_edges: #This script returns nodes and edges that pass the desired level of correlation and pvalue. It also calculates modularity, and the Zi-Pi calculations
   version: "1.0"
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      pvalues=rules.pvalues.output,
      corr=rules.fastspar.output
   output:
      nodes=report("data/network/{sample}–{group}/nodes–{groupcategory}.tsv",
      caption="../report/network_nodes_edges.rst",
      category="SparCC",
      subcategory="Nodes and edges",
      labels={
         "File type": "Nodes - Text file",
         "Pairwise comparison": "{groupcategory}",
         "Grouping category": "{group}"}),
      edges=report("data/network/{sample}–{group}/edges–{groupcategory}.tsv",
      caption="../report/network_nodes_edges.rst",
      category="SparCC",
      subcategory="Nodes and edges",
      labels={
         "File type": "Edges - Text file",
         "Pairwise comparison": "{groupcategory}",
         "Grouping category": "{group}"})
   params:
      pvalue=config["sparcc_pvalue"][0],
      threshold=config["sparcc_corr"][0]
   message: "Creating nodes and edges for the calculated SparCC for {wildcards.sample} {wildcards.group} {wildcards.groupcategory}"
   shell:
      "Rscript --vanilla ./workflow/scripts/nodes_edges.R -i {input.corr} -p {input.pvalues} -n {output.nodes} -e {output.edges} -a {params.threshold} -b {params.pvalue} >/dev/null"

rule zipi: # This script plots the Zi-Pi plots from the nodes_and_edges output
   version: "1.0"
   conda:
      "../../workflow/envs/ggrepel.yaml"
   input:
      rules.nodes_and_edges.output.nodes
   output:
      report("data/plots/network/ZiPi_{sample}–{group}/ZiPi-{groupcategory}.svg",
      caption="../report/network_zipi.rst",
      category="SparCC",
      subcategory="ZiPi plots",
      labels={
         "File type": "ZiPi plot",
         "Comparison": "{groupcategory}",
         "Grouping category": "{group}"})
   params:
      pvalue=config["sparcc_pvalue"][0],
      corr=config["sparcc_corr"][0],
      core=config["threshold"][0]
   message: "Generating Zi-Pi plots for {wildcards.sample} {wildcards.group} {wildcards.groupcategory}"
   shell:
      "Rscript --vanilla ./workflow/scripts/zipi_graph.R -i {input} -o {output} -t {params.corr} -c {params.core} -p {params.pvalue}  >/dev/null"


def zipi_output_aggregation(wildcards):
    checkpoint_output = checkpoints.makeCore.get(**wildcards).output[0]    
    return expand("data/plots/network/ZiPi_{sample}–{group}/ZiPi-{groupcategory}.svg",
           sample=wildcards.sample,
           group=wildcards.group,
           groupcategory=glob_wildcards(os.path.join(checkpoint_output, "{groupcategory}.tsv")).groupcategory)


rule network_cleanup: # Cleans up the temporary files from the network calculations
   version: "1.0"
   input:
      zipi_output_aggregation
   output:
      temporary(touch("tmp/network_{sample}–{group}.done"))
   message: "Cleaning up after networks of {wildcards.sample} variable {wildcards.group}"
   shell:
      "rm -rf data/network/{wildcards.sample}–{wildcards.group}_corr && "
      "rm -rf data/network/{wildcards.sample}–{wildcards.group}/bootstrap* "
      
rule network:
   input:
      expand("tmp/network_{sample}–{group}.done",  sample=config["mysample"],group=config["group"])