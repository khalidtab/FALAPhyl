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
   shell:
      "mkdir -p data/network/ &&"
      "csvcut {input} -t -c {params.group} | uniq | tail -n +2 > {output}"

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
   shell:
      "mkdir -p {output} &&"
      "Rscript --vanilla ./workflow/scripts/make_core.R -i {input.biom} -m {input.map} -c {params.group} -t {params.threshold} -o {output}"

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
      "fastspar -c {input.coreFolder}core+$i.tsv -r {output.sparcc}corr+$i.tsv -a {output.cov}cov+$i.tsv -e {params.corr}"
      "; done' > temp/fastspar_on_core_{wildcards.sample}.sh && "
      "chmod +x temp/fastspar_on_core_{wildcards.sample}.sh &&"
      "bash temp/fastspar_on_core_{wildcards.sample}.sh >/dev/null"

rule bootstrap:
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      coreFolder=rules.makeCore.output,
      mycat="data/network/{sample}_categories.txt",
      sparcc=rules.fastspar_on_core.output.sparcc
   output:
      bootstrapDone=temporary(touch("temp/bootstrap_{sample}.done"))
   params:
      bootstrap=config["sparcc_bootstrap"][0]
   shell:
      "echo 'for i in $(cat {input.mycat}); do  "
      "mkdir -p data/network/{wildcards.sample}/bootstrap+$i/fake_data/ && "
      "fastspar_bootstrap --otu_table {input.coreFolder}core+$i.tsv --number {params.bootstrap} --prefix data/network/{wildcards.sample}/bootstrap+$i/fake_data/"
      "; done' > temp/bootstrap_{wildcards.sample}.sh && "
      "chmod +x temp/bootstrap_{wildcards.sample}.sh &&"
      "bash temp/bootstrap_{wildcards.sample}.sh >/dev/null"
  

rule bootstrap_fastspar:
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      coreFolder=rules.bootstrap.output.bootstrapDone,
      mycat="data/network/{sample}_categories.txt"
   output:
      fastspar=temporary(touch("temp/fastspar_{sample}.done"))
   params:
      bootstrap=config["sparcc_bootstrap"][0]
   shell:
      "echo 'for i in $(cat {input.mycat}); do  "
      "mkdir -p data/network/{wildcards.sample}/bootstrap+$i && "
      "for y in {{0..{params.bootstrap}}}; do "
      "if [ $y -lt {params.bootstrap} ]; then fastspar --otu_table data/network/{wildcards.sample}/bootstrap+$i/fake_data/_$y.tsv -i 5 --correlation data/network/{wildcards.sample}/bootstrap+$i/fake_data_corr_$y.tsv --covariance data/network/{wildcards.sample}/bootstrap+$i/fake_data_cov_$y.tsv"
      "; fi ; done; done' > temp/fastspar_bootstrap_{wildcards.sample}.sh && "
      "chmod +x temp/fastspar_bootstrap_{wildcards.sample}.sh && "
      "bash temp/fastspar_bootstrap_{wildcards.sample}.sh >/dev/null"
  


rule pvalues:
   version: "1.0"
   conda:
      "../../workflow/envs/fastspar.yaml"
   input:
      coreFolder=rules.bootstrap_fastspar.output,
      mycat="data/network/{sample}_categories.txt",
      folder=rules.makeCore.output.network
   output:
      fastspar=temporary(touch("temp/pvalues_{sample}.done"))
   params:
      bootstrap=config["sparcc_bootstrap"][0]
   shell:
      "echo 'for i in $(cat {input.mycat}); do "
      "fastspar_pvalues --otu_table data/network/{wildcards.sample}/core+$i.tsv --correlation data/network/{wildcards.sample}_corr/corr+$i.tsv --prefix data/network/{wildcards.sample}/bootstrap+$i/fake_data_corr_ --permutations {params.bootstrap} --outfile data/network/{wildcards.sample}_corr/pvalue+$i.tsv "
      "; done' > temp/fastspar_{wildcards.sample}.sh && "
      "chmod +x temp/fastspar_{wildcards.sample}.sh && "
      "bash temp/fastspar_{wildcards.sample}.sh >/dev/null"

rule nodes_and_edges:
   version: "1.0"
   conda:
      "../../workflow/envs/R_with_graphing.yaml"
   input:
      pvalues=rules.pvalues.output,
      mycat="data/network/{sample}_categories.txt"
   output:
      temporary(touch("temp/table_{sample}.done"))
   params:
      pvalue=config["sparcc_pvalue"][0],
      threshold=config["threshold"][0]
   shell:
      "mkdir -p data/network/{wildcards.sample} && echo 'for i in $(cat {input.mycat}); do "
      "Rscript --vanilla ./workflow/scripts/gephi_zipi.R -i data/network/{wildcards.sample}_corr/corr+$i.tsv -p data/network/{wildcards.sample}_corr/pvalue+$i.tsv -n data/network/{wildcards.sample}/nodes+$i.tsv -e data/network/{wildcards.sample}/edges+$i.tsv -a {params.threshold} -b {params.pvalue}"
      "; done' > temp/network_{wildcards.sample}.sh && "
      "chmod +x temp/network_{wildcards.sample}.sh && "
      "bash temp/network_{wildcards.sample}.sh"


rule run_fastspar:
   version: "1.0"
   conda:
      "../../workflow/envs/NMDS.yaml"
   input:
      pvalues=rules.nodes_and_edges.output,
      mycat="data/network/{sample}_categories.txt"
   output:
      temporary(touch("temp/zipiGraph_{sample}.done"))
   shell:
      "mkdir -p data/logs data/plots && echo 'for i in $(cat {input.mycat}); do "
      "Rscript --vanilla ./workflow/scripts/zipi_graph.R -i data/network/{wildcards.sample}/nodes+$i.tsv -o data/plots/{wildcards.sample}_zipi+$i.json &&"
      "xvfb-run --auto-servernum orca graph data/plots/{wildcards.sample}_zipi+$i.json -o data/plots/zipi_{wildcards.sample}+$i.svg -f svg 2>> data/logs/zipi_{wildcards.sample}+$x.log || true &&"
      "rm data/plots/{wildcards.sample}_zipi+$i.json &&"
      "rm -rf data/network/{wildcards.sample}_corr"
      "; done' > temp/zipi_{wildcards.sample}.sh &&"
      "chmod +x temp/zipi_{wildcards.sample}.sh && "
      "bash temp/zipi_{wildcards.sample}.sh"

