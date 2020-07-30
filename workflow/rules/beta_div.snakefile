rule beta_div:
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      "data/biom/{sample}.qza"
   params: 
      dist=expand("{dist}",dist=config["distances"])
   output:
      qza=temporary(expand("data/distance/beta_div/{{sample}}+{dist}.qza",dist=config["distances"])),
      qza_pcoa=temporary(expand("data/distance/PCoA/PCoA_{{sample}}+{dist}.qza",dist=config["distances"])),
      tsv=expand("data/distance/beta_div/{{sample}}+{dist}.tsv",dist=config["distances"]),
      tsv_pcoa=expand("data/distance/PCoA/PCoA_{{sample}}+{dist}.tsv",dist=config["distances"])
   shell: 
      "mkdir -p data/distance/beta_div data/distance/PCoA &&"
      "echo 'for x in {params.dist}; do qiime diversity beta --i-table data/biom/{wildcards.sample}.qza --o-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --p-metric $x &&"
      "qiime diversity pcoa --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --o-pcoa data/distance/PCoA/PCoA_{wildcards.sample}+$x.qza &&"
      "qiime tools export --input-path data/distance/beta_div/{wildcards.sample}+$x.qza --output-path data/distance/beta_div/{wildcards.sample}+$x &&"
      "qiime tools export --input-path data/distance/PCoA/PCoA_{wildcards.sample}+$x.qza --output-path data/distance/PCoA/PCoA_{wildcards.sample}+$x &&"
      "cp data/distance/beta_div/{wildcards.sample}+$x/distance-matrix.tsv data/distance/beta_div/{wildcards.sample}+$x.tsv &&"
      "cp data/distance/PCoA/PCoA_{wildcards.sample}+$x/ordination.txt data/distance/PCoA/PCoA_{wildcards.sample}+$x.tsv && "
      "rm -rf data/distance/PCoA/PCoA_{wildcards.sample}+$x/ && rm -rf data/distance/beta_div/{wildcards.sample}+$x/"
      "; done' > tmp/beta_div_PCoA_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_PCoA_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_PCoA_{wildcards.sample}.sh"


rule nmds:
   version: "1.0"
   conda: "../../workflow/envs/NMDS.yaml"
   input:
      rules.beta_div.output.tsv
   output: 
      temporary(expand("data/plots/NMDS_{{sample}}+{dist}.json",dist=config["distances"])),
      report(expand("data/plots/NMDS_{{sample}}+{dist}_1.svg",dist=config["distances"]))
   params:
      dist=expand("{dist}",dist=config["distances"]),
      group=config["group"][0],
      color=config["color"][0]
   log:
      expand("data/logs/NMDS_{{sample}}+{dist}.log", dist=config["distances"])
   shell:
      "echo 'for x in {params.dist}; do Rscript --vanilla ./workflow/scripts/NMDS.R -i data/distance/beta_div/{wildcards.sample}+$x.tsv -o data/plots/NMDS_{wildcards.sample}+$x.json -m data/map/{wildcards.sample}.txt -g {params.group} -c {params.color} > data/logs/NMDS_{wildcards.sample}+$x.log && "
      "xvfb-run --auto-servernum orca graph data/plots/NMDS_{wildcards.sample}+$x.json -o data/plots/NMDS_{wildcards.sample}+$x.svg -f svg 2>> data/logs/NMDS_{wildcards.sample}+$x.log || true; done' > tmp/beta_div_NMDS_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_NMDS_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_NMDS_{wildcards.sample}.sh"

rule nmds_hull:
   version: "1.0"
   conda: "../../workflow/envs/NMDS.yaml"
   input:
      rules.beta_div.output.tsv
   output: 
      report(expand("data/plots/NMDS_hull_{{sample}}+{dist}_1.svg",dist=config["distances"])),
      temporary(expand("data/plots/NMDS_hull_{{sample}}+{dist}.json",dist=config["distances"]))
   params:
      dist=expand("{dist}",dist=config["distances"]),
      group=config["group"][0],
      color=config["color"][0]
   log:
      expand("data/logs/NMDS_hull_{{sample}}+{dist}.log", dist=config["distances"])
   shell:
      "echo 'for x in {params.dist}; do Rscript --vanilla ./workflow/scripts/NMDS_hull.R -i data/distance/beta_div/{wildcards.sample}+$x.tsv -o data/plots/NMDS_hull_{wildcards.sample}+$x.json -m data/map/{wildcards.sample}.txt -g {params.group} -c {params.color} > data/logs/NMDS_hull_{wildcards.sample}+$x.log && xvfb-run --auto-servernum orca graph data/plots/NMDS_hull_{wildcards.sample}+$x.json -o data/plots/NMDS_hull_{wildcards.sample}+$x.svg -f svg 2>> data/logs/NMDS_hull_{wildcards.sample}+$x.log || true; done' > tmp/beta_div_NMDS_hull_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_NMDS_hull_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_NMDS_hull_{wildcards.sample}.sh"

rule get_colors:
   version: "1.0"
   conda:
      "../../workflow/envs/csvkit.yaml"
   input:
      "data/map/{sample}.txt"
   params:
      color=config["color"][0]
   output:
      temporary("data/map/{sample}_color.txt")
   shell:
      "csvcut {input} -t -c {params.color} | uniq | tail -n +2 > {output}"

rule pcoa_svg:
   version: "1.0"
   conda:
      "../../workflow/envs/phylotoast.yaml"
   input:
      betaDiv=rules.beta_div.output.tsv_pcoa,
      color=rules.get_colors.output
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=config["group"][0]
   output:
      report(expand("data/plots/PCoA_{{sample}}+{dist}.svg",dist=config["distances"]))
   shell: 
      "mkdir -p data/plots &&"
      "echo 'for x in {params.dist}; do xvfb-run --auto-servernum PCoA.py -i data/distance/PCoA/PCoA_{wildcards.sample}+$x.tsv -m data/map/{wildcards.sample}.txt -b {params.group} -d 2 -c {input.color} -o data/plots/PCoA_{wildcards.sample}+$x.svg ; done'"
      "> tmp/SVG_PCoA_{wildcards.sample}.sh &&"
      "chmod +x tmp/SVG_PCoA_{wildcards.sample}.sh &&"
      "bash tmp/SVG_PCoA_{wildcards.sample}.sh "