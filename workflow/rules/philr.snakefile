#rule correlation_clustering:
#   version: "1.0"
#   conda:
#      "../../workflow/envs/qiime2.yaml"
#   input:
#      "data/biom/{sample}.biom"
#   output:
#      "data/tree/{sample}.tre"
#   message: "Generating correlation tree for {wildcards.sample}"
#   shell:
#      "mkdir -p data/tsv &&"
#      "biom convert -i {input} -o data/tsv/{wildcards.sample}_for_tree.tsv --to-tsv --header-key taxonomy &&"
#      'biom convert -i data/tsv/{wildcards.sample}_for_tree.tsv -o data/biom/{wildcards.sample}_for_tree.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy &&'
#      "qiime tools import --input-path data/biom/{wildcards.sample}_for_tree.biom --output-path data/biom/{wildcards.sample}_for_tree.qza --type FeatureTable[Frequency] --input-format BIOMV100Format >/dev/null &&"
#      "qiime gneiss correlation-clustering --i-table data/biom/{wildcards.sample}_for_tree.qza --o-clustering data/tree/{wildcards.sample}.qza  >/dev/null &&"
#      "qiime tools export --input-path data/tree/{wildcards.sample}.qza --output-path data/tree/tree_{wildcards.sample}  >/dev/null &&"
#      "mv data/tree/tree_{wildcards.sample}/tree.nwk {output} &&"
#      "rm data/tree/{wildcards.sample}.qza data/biom/{wildcards.sample}_for_tree.qza data/tsv/{wildcards.sample}_for_tree.tsv &&"
#      "rm -rf data/tree/tree_{wildcards.sample}"




rule philr_distances:
   version: "1.0"
   conda:
      "../../workflow/envs/philr.yaml"
   input:
      biom = "data/biom/{sample}.biom",
      tree = "data/tree/{sample}.tre",
      map = "data/map/{sample}.txt"
   params: 
      group=config["group"][0]
   output: "data/distance/beta_div/{sample}+philr.tsv"
   message: "Calculating PhILR distances for {wildcards.sample}"
   shell:
      "mkdir -p data/tsv &&"
      "Rscript --vanilla ./workflow/scripts/philr.R -i {input.biom} -t {input.tree} -m {input.map} -g {params.group} -o {output} 2>> /dev/null"
      
rule philr_nmds:
   version: "1.0"
   conda: "../../workflow/envs/NMDS.yaml"
   input:
      rules.philr_distances.output
   output: 
      temporary("data/plots/NMDS_{{sample}}+philr.json")
   params:
      group=config["group"][0],
      color=config["color"][0]
   log:
      "data/logs/NMDS_{{sample}}+philr.log"
   message: "Generating NMDS plots using PhILR distances for {wildcards.sample}"
   shell:
      "Rscript --vanilla ./workflow/scripts/NMDS.R -i data/distance/beta_div/{wildcards.sample}+PhILR.tsv -o data/plots/NMDS_{wildcards.sample}+PhILR.json -m data/map/{wildcards.sample}.txt -g {params.group} -c {params.color} > data/logs/NMDS_{wildcards.sample}+PhILR.log 2>>/dev/null && "
      "xvfb-run --auto-servernum orca graph data/plots/NMDS_{wildcards.sample}+PhILR.json -o data/plots/NMDS_{wildcards.sample}+PhILR.svg -f svg 2>> data/logs/NMDS_{wildcards.sample}+PhILR.log || true > tmp/beta_div_NMDS_Philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_NMDS_Philr_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_NMDS_Philr_{wildcards.sample}.sh"