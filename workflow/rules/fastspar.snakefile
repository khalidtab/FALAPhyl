rule get_categories: # Prepares the chosen categories
   version: "1.0"
   conda:
      "../../workflow/envs/csvkit.yaml"
   input:
      "data/map/{sample}.txt"
   params:
      group=expand("{group}",group=config["group"])
   output:
      temporary("data/network/{sample}_categories.txt")
   message: "Extracting catergories from the mapping file for {wildcards.sample}"
   shell:
      "mkdir -p data/network/ &&"
      "echo 'for x in {params.group}; do "
      "csvcut {input} -t -c $x | tail -n +2 | sort -u | csvformat -T >> data/network/{wildcards.sample}_categories.txt; done ' > tmp/getCategories_{wildcards.sample}.sh &&"
      "chmod +x tmp/getCategories_{wildcards.sample}.sh &&"
      "bash tmp/getCategories_{wildcards.sample}.sh"


rule makeCore: # SparCC works better when the feature table has less zeros. Therefore, this script will filter the feature tables to the desired level of "coreness" based on the input file values
   version: "1.0"
   conda:
      "../../workflow/envs/phyloseq.yaml"
   input:
      biom="data/biom/{sample}_temp.biom",
      map="data/map/{sample}.txt"
   params:
      threshold=config["threshold"][0],
      group=config["group"]
   output:
      network=directory("data/network/{sample}/")
   message: "Making core files for {wildcards.sample}"
   shell:
      "mkdir -p {output} data/logs &&"
      "echo 'for x in {params.group}; do "
      "Rscript --vanilla ./workflow/scripts/make_core.R -i {input.biom} -m {input.map} -c $x -t {params.threshold} -o {output} ; done ' > tmp/makeCores_{wildcards.sample}.sh &&"
      "chmod +x tmp/makeCores_{wildcards.sample}.sh &&"
      "bash tmp/makeCores_{wildcards.sample}.sh"


rule fastspar_on_core: # SparCC (through the fastpar algorithm) will be done on the core feature table
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      coreFolder=rules.makeCore.output,
      mycat="data/network/{sample}_categories.txt"
   output:
      sparcc=directory("data/network/{sample}_corr/"),
      cov=temporary(directory("data/network/{sample}_cov/"))
   params:
      corr=config["sparcc_corr"][0]
   message: "Calculating SparCC on core of {wildcards.sample}"
   shell:
      "mkdir -p {output.sparcc} {output.cov} &&"
      "echo 'for i in $(cat {input.mycat}); do "
      "fastspar -c {input.coreFolder}_core+$i.tsv -r {output.sparcc}/corr+$i.tsv -a {output.cov}/cov+$i.tsv -e {params.corr} >/dev/null"
      "; done' > tmp/fastspar_on_core_{wildcards.sample}.sh && "
      "chmod +x tmp/fastspar_on_core_{wildcards.sample}.sh &&"
      "bash tmp/fastspar_on_core_{wildcards.sample}.sh"

rule bootstrap: # To calculate P-values, bootstraps need to be done from the core feature table. This script accomplishes this.
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      coreFolder=rules.makeCore.output,
      mycat="data/network/{sample}_categories.txt",
      sparcc=rules.fastspar_on_core.output.sparcc
   output:
      bootstrapDone=temporary(touch("tmp/bootstrap_{sample}.done"))
   params:
      bootstrap=config["sparcc_bootstrap"][0]
   message: "Making bootstrap on {wildcards.sample}"
   shell:
      "echo 'for i in $(cat {input.mycat}); do  "
      "mkdir -p data/network/{wildcards.sample}/bootstrap+$i/fake_data/ && "
      "fastspar_bootstrap --otu_table {input.coreFolder}_core+$i.tsv --number {params.bootstrap} --prefix data/network/{wildcards.sample}/bootstrap+$i/fake_data/  >/dev/null"
      "; done' > tmp/bootstrap_{wildcards.sample}.sh && "
      "chmod +x tmp/bootstrap_{wildcards.sample}.sh &&"
      "bash tmp/bootstrap_{wildcards.sample}.sh"
  

rule bootstrap_fastspar: # TO calculate p-values, correlations on the bootstraps need to be done. This script accomplishes this.
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      coreFolder=rules.bootstrap.output.bootstrapDone,
      mycat="data/network/{sample}_categories.txt"
   output:
      fastspar=temporary(touch("tmp/fastspar_{sample}.done"))
   params:
      bootstrap=config["sparcc_bootstrap"][0]
   message: "Calculating SparCC on bootstraps of {wildcards.sample}"
   shell:
      "echo 'for i in $(cat {input.mycat}); do  "
      "mkdir -p data/network/{wildcards.sample}/bootstrap+$i && "
      "for y in {{0..{params.bootstrap}}}; do "
      "if [ $y -lt {params.bootstrap} ]; then fastspar --otu_table data/network/{wildcards.sample}/bootstrap+$i/fake_data/_$y.tsv -i 5 --correlation data/network/{wildcards.sample}/bootstrap+$i/fake_data_corr_$y.tsv --covariance data/network/{wildcards.sample}/bootstrap+$i/fake_data_cov_$y.tsv  >/dev/null"
      "; fi ; done; done' > tmp/fastspar_bootstrap_{wildcards.sample}.sh && "
      "chmod +x tmp/fastspar_bootstrap_{wildcards.sample}.sh && "
      "bash tmp/fastspar_bootstrap_{wildcards.sample}.sh"

rule pvalues: # Calculates p-values by using the correlations from the core feature table, and the boostraps
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      coreFolder=rules.bootstrap_fastspar.output,
      mycat="data/network/{sample}_categories.txt",
      folder=rules.makeCore.output.network
   output:
      fastspar=temporary(touch("tmp/pvalues_{sample}.done"))
   params:
      bootstrap=config["sparcc_bootstrap"][0]
   message: "pvalues calculation for the features in networks of {wildcards.sample}"
   shell:
      "echo 'for i in $(cat {input.mycat}); do "
      "fastspar_pvalues --otu_table data/network/{wildcards.sample}_core+$i.tsv --correlation data/network/{wildcards.sample}_corr/corr+$i.tsv --prefix data/network/{wildcards.sample}/bootstrap+$i/fake_data_corr_ --permutations {params.bootstrap} --outfile data/network/{wildcards.sample}_corr/pvalue+$i.tsv  >/dev/null"
      "; done' > tmp/fastspar_{wildcards.sample}.sh && "
      "chmod +x tmp/fastspar_{wildcards.sample}.sh && "
      "bash tmp/fastspar_{wildcards.sample}.sh"

rule nodes_and_edges: #This script returns nodes and edges that pass the desired level of correlation and pvalue. It also calculates modularity, and the Zi-Pi calculations
   version: "1.0"
   conda:
      "../../workflow/envs/phyloseq.yaml"
   input:
      pvalues=rules.pvalues.output,
      mycat="data/network/{sample}_categories.txt"
   output:
      temporary(touch("tmp/table_{sample}.done"))
   params:
      pvalue=config["sparcc_pvalue"][0],
      threshold=config["sparcc_corr"][0]
   message: "Creating nodes and edges for the calculated SparCC for {wildcards.sample}"
   shell:
      "mkdir -p data/network/{wildcards.sample} && echo 'for i in $(cat {input.mycat}); do "
      "Rscript --vanilla ./workflow/scripts/nodes_edges.R -i data/network/{wildcards.sample}_corr/corr+$i.tsv -p data/network/{wildcards.sample}_corr/pvalue+$i.tsv -n data/network/{wildcards.sample}/nodes+$i.tsv -e data/network/{wildcards.sample}/edges+$i.tsv -a {params.threshold} -b {params.pvalue} >/dev/null"
      "; done' > tmp/network_{wildcards.sample}.sh && "
      "chmod +x tmp/network_{wildcards.sample}.sh && "
      "bash tmp/network_{wildcards.sample}.sh"

rule zipi: # This script plots the Zi-Pi plots from the nodes_and_edges output
   version: "1.0"
   conda:
      "../../workflow/envs/ggrepel.yaml"
   input:
      pvalues=rules.nodes_and_edges.output,
      mycat="data/network/{sample}_categories.txt"
   output:
      temporary(touch("tmp/zipiGraph_{sample}.done"))
   params:
      pvalue=config["sparcc_pvalue"][0],
      corr=config["sparcc_corr"][0],
      core=config["threshold"][0]
   message: "Generating Zi-Pi plots for {wildcards.sample}"
   shell:
      "mkdir -p data/logs data/plots && echo 'for i in $(cat {input.mycat}); do "
      "Rscript --vanilla ./workflow/scripts/zipi_graph.R -i data/network/{wildcards.sample}/nodes+$i.tsv -o data/plots/ZiPi_{wildcards.sample}+$i.svg -t {params.corr} -c {params.core} -p {params.pvalue}  >/dev/null"
      "; done' > tmp/zipi_{wildcards.sample}.sh &&"
      "chmod +x tmp/zipi_{wildcards.sample}.sh && "
      "bash tmp/zipi_{wildcards.sample}.sh >/dev/null"

rule network: # Cleans up the temporary files from the network calculations
   version: "1.0"
   input:
      zipi=rules.zipi.output,
      nodes=rules.nodes_and_edges.output
   output:
      temporary(touch("tmp/network_{sample}.done"))
   message: "Cleaning up after networks of {wildcards.sample}"
   shell:
      "rm -rf data/network/{wildcards.sample}_corr &&"
      "rm -rf data/network/{wildcards.sample}/bootstrap*"