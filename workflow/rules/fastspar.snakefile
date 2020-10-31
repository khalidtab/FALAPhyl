rule get_categories:
   version: "1.0"
   conda:
      "../../workflow/envs/csvkit.yaml"
   input:
      "data/map/{sample}.txt"
   params:
      group=config["group"][0]
   output:
      temporary("data/network/{sample}_categories.txt")
   message: "Extracting catergories from the mapping file for {wildcards.sample}"
   shell:
      "mkdir -p data/network/ &&"
      "csvcut {input} -t -c {params.group} | sort -u | head -n -1 > {output}"

rule makeCore:
   version: "1.0"
   conda:
      "../../workflow/envs/R_with_graphing.yaml"
   input:
      biom="data/biom/{sample}_temp.biom",
      map="data/map/{sample}.txt"
   params:
      threshold=config["threshold"][0],
      group=config["group"][0]
   output:
      network=temporary(directory("data/network/{sample}/"))
   message: "Making core files for {wildcards.sample}"
   shell:
      "mkdir -p {output} &&"
      "Rscript --vanilla ./workflow/scripts/make_core.R -i {input.biom} -m {input.map} -c {params.group} -t {params.threshold} -o {output}  >/dev/null"

rule fastspar_on_core:
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
      "fastspar -c {input.coreFolder}core+$i.tsv -r {output.sparcc}corr+$i.tsv -a {output.cov}cov+$i.tsv -e {params.corr} >/dev/null"
      "; done' > tmp/fastspar_on_core_{wildcards.sample}.sh && "
      "chmod +x tmp/fastspar_on_core_{wildcards.sample}.sh &&"
      "bash tmp/fastspar_on_core_{wildcards.sample}.sh"

rule bootstrap:
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
      "fastspar_bootstrap --otu_table {input.coreFolder}core+$i.tsv --number {params.bootstrap} --prefix data/network/{wildcards.sample}/bootstrap+$i/fake_data/  >/dev/null"
      "; done' > tmp/bootstrap_{wildcards.sample}.sh && "
      "chmod +x tmp/bootstrap_{wildcards.sample}.sh &&"
      "bash tmp/bootstrap_{wildcards.sample}.sh"
  

rule bootstrap_fastspar:
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

rule pvalues:
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
   message: "pvalues calculation for the features of {wildcards.sample}"
   shell:
      "echo 'for i in $(cat {input.mycat}); do "
      "fastspar_pvalues --otu_table data/network/{wildcards.sample}/core+$i.tsv --correlation data/network/{wildcards.sample}_corr/corr+$i.tsv --prefix data/network/{wildcards.sample}/bootstrap+$i/fake_data_corr_ --permutations {params.bootstrap} --outfile data/network/{wildcards.sample}_corr/pvalue+$i.tsv  >/dev/null"
      "; done' > tmp/fastspar_{wildcards.sample}.sh && "
      "chmod +x tmp/fastspar_{wildcards.sample}.sh && "
      "bash tmp/fastspar_{wildcards.sample}.sh"

rule nodes_and_edges:
   version: "1.0"
   conda:
      "../../workflow/envs/R_with_graphing.yaml"
   input:
      pvalues=rules.pvalues.output,
      mycat="data/network/{sample}_categories.txt"
   output:
      temporary(touch("tmp/table_{sample}.done"))
   params:
      pvalue=config["sparcc_pvalue"][0],
      threshold=config["threshold"][0]
   message: "Creating nodes and edges for the calculated SparCC for {wildcards.sample}"
   shell:
      "mkdir -p data/network/{wildcards.sample} && echo 'for i in $(cat {input.mycat}); do "
      "Rscript --vanilla ./workflow/scripts/gephi_zipi.R -i data/network/{wildcards.sample}_corr/corr+$i.tsv -p data/network/{wildcards.sample}_corr/pvalue+$i.tsv -n data/network/{wildcards.sample}/nodes+$i.tsv -e data/network/{wildcards.sample}/edges+$i.tsv -a {params.threshold} -b {params.pvalue} >/dev/null"
      "; done' > tmp/network_{wildcards.sample}.sh && "
      "chmod +x tmp/network_{wildcards.sample}.sh && "
      "bash tmp/network_{wildcards.sample}.sh"

rule network:
   version: "1.0"
   conda:
      "../../workflow/envs/NMDS.yaml"
   input:
      pvalues=rules.nodes_and_edges.output,
      mycat="data/network/{sample}_categories.txt"
   output:
      temporary(touch("tmp/zipiGraph_{sample}.done"))
   message: "Generating Zi-Pi plots for {wildcards.sample}"
   shell:
      "mkdir -p data/logs data/plots && echo 'for i in $(cat {input.mycat}); do "
      "Rscript --vanilla ./workflow/scripts/zipi_graph.R -i data/network/{wildcards.sample}/nodes+$i.tsv -o data/plots/{wildcards.sample}_zipi+$i.json  >/dev/null &&"
      "xvfb-run --auto-servernum orca graph data/plots/{wildcards.sample}_zipi+$i.json -o data/plots/zipi_{wildcards.sample}+$i.svg -f svg 2>> data/logs/zipi_{wildcards.sample}+$x.log || true &&"
      "rm data/plots/{wildcards.sample}_zipi+$i.json &&"
      "rm -rf data/network/{wildcards.sample}_corr"
      "; done' > tmp/zipi_{wildcards.sample}.sh &&"
      "chmod +x tmp/zipi_{wildcards.sample}.sh && "
      "bash tmp/zipi_{wildcards.sample}.sh >/dev/null"