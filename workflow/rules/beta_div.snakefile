rule get_colors:
   version: "1.0"
   conda:
      "../../workflow/envs/csvkit.yaml"
   input:
      "data/map/{sample}.txt"
   params:
      group=expand("{group}",group=config["group"]),
      color=config["color"][0]
   output:
      temporary(expand("data/map/color_{{sample}}+{group}.txt",group=config["group"]))
   message: "Extracting colors from the mapping file for {wildcards.sample}"
   shell:
      "echo 'for x in {params.group}; do "
      "y=$(printf \"%s{params.color}\" $x) &&" 
      "csvcut {input} -t -c $x,$y | tail -n +2 | sort -u | csvformat -T | csvcut -t -c 2 > data/map/color_{wildcards.sample}+$x.txt; done ' > tmp/getColors_{wildcards.sample}.sh &&"
      "chmod +x tmp/getColors_{wildcards.sample}.sh &&"
      "bash tmp/getColors_{wildcards.sample}.sh"

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
   message: "Calculating beta diversity for {wildcards.sample}"
   shell: 
      "mkdir -p data/distance/beta_div data/distance/PCoA &&"
      "echo 'for x in {params.dist}; do qiime diversity beta --i-table data/biom/{wildcards.sample}.qza --o-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --p-metric $x  >/dev/null &&"
      "qiime diversity pcoa --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --o-pcoa data/distance/PCoA/PCoA_{wildcards.sample}+$x.qza  >/dev/null &&"
      "qiime tools export --input-path data/distance/beta_div/{wildcards.sample}+$x.qza --output-path data/distance/beta_div/{wildcards.sample}+$x  >/dev/null &&"
      "qiime tools export --input-path data/distance/PCoA/PCoA_{wildcards.sample}+$x.qza --output-path data/distance/PCoA/PCoA_{wildcards.sample}+$x  >/dev/null &&"
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
      temporary(expand("data/plots/NMDS_{{sample}}+{dist}+{group}.json",dist=config["distances"],group=config["group"])),
      report(expand("data/plots/NMDS_{{sample}}+{dist}+{group}_1.svg",dist=config["distances"],group=config["group"]))
   params:
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"]),
      color=config["color"][0]
   log:
      expand("data/logs/NMDS_{{sample}}+{dist}+{group}.log", dist=config["distances"],group=config["group"])
   message: "Generating NMDS plots for {wildcards.sample}"
   shell:
      "echo 'for x in {params.dist}; do for w in {params.group}; do "
      "y=$(printf \"%s{params.color}\" $w) && "
      " Rscript --vanilla ./workflow/scripts/NMDS.R -i data/distance/beta_div/{wildcards.sample}+$x.tsv -o data/plots/NMDS_{wildcards.sample}+$x+$w.json -m data/map/{wildcards.sample}.txt -g $w -c $y > data/logs/NMDS_{wildcards.sample}+$x+$w.log 2>>/dev/null && "
      "xvfb-run --auto-servernum orca graph data/plots/NMDS_{wildcards.sample}+$x+$w.json -o data/plots/NMDS_{wildcards.sample}+$x+$w.svg -f svg 2>> data/logs/NMDS_{wildcards.sample}+$x+$w.log || true; done; done' > tmp/beta_div_NMDS_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_NMDS_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_NMDS_{wildcards.sample}.sh"

rule nmds_hull:
   version: "1.0"
   conda: "../../workflow/envs/NMDS.yaml"
   input:
      rules.beta_div.output.tsv
   output: 
      temporary(expand("data/plots/NMDS_HULL_{{sample}}+{dist}+{group}.json",dist=config["distances"],group=config["group"])),
      report(expand("data/plots/NMDS_HULL_{{sample}}+{dist}+{group}_1.svg",dist=config["distances"],group=config["group"]))
   params:
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"]),
      color=config["color"][0]
   log:
      expand("data/logs/NMDS_HULL_{{sample}}+{dist}+{group}.log", dist=config["distances"],group=config["group"])
   message: "Generating NMDS_HULL plots for {wildcards.sample}"
   shell:
      "echo 'for x in {params.dist}; do for w in {params.group}; do "
      "y=$(printf \"%s{params.color}\" $w) && "
      " Rscript --vanilla ./workflow/scripts/NMDS_hull.R -i data/distance/beta_div/{wildcards.sample}+$x.tsv -o data/plots/NMDS_HULL_{wildcards.sample}+$x+$w.json -m data/map/{wildcards.sample}.txt -g $w -c $y > data/logs/NMDS_HULL_{wildcards.sample}+$x+$w.log 2>>/dev/null && "
      "xvfb-run --auto-servernum orca graph data/plots/NMDS_HULL_{wildcards.sample}+$x+$w.json -o data/plots/NMDS_HULL_{wildcards.sample}+$x+$w.svg -f svg 2>> data/logs/NMDS_HULL_{wildcards.sample}+$x+$w.log || true; done; done' > tmp/beta_div_NMDS_HULL_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_NMDS_HULL_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_NMDS_HULL_{wildcards.sample}.sh"


rule pcoa_svg:
   version: "1.0"
   conda:
      "../../workflow/envs/phylotoast.yaml"
   input:
      betaDiv=rules.beta_div.output.tsv_pcoa,
      color=rules.get_colors.output
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"])
   output:
      report(expand("data/plots/PCoA_{{sample}}+{dist}+{group}.svg",dist=config["distances"],group=config["group"]))
   message: "Generating PCoA plots for {wildcards.sample}"
   shell: 
      "mkdir -p data/plots &&"
      "echo 'for x in {params.dist}; do for y in {params.group};" 
      "do xvfb-run --auto-servernum PCoA.py -i data/distance/PCoA/PCoA_{wildcards.sample}+$x.tsv -m data/map/{wildcards.sample}.txt -b $y -d 2 -c data/map/color_{wildcards.sample}+$y.txt -o data/plots/PCoA_{wildcards.sample}+$x+$y.svg ; done ; done'"
      "> tmp/SVG_PCoA_{wildcards.sample}.sh &&"
      "chmod +x tmp/SVG_PCoA_{wildcards.sample}.sh &&"
      "bash tmp/SVG_PCoA_{wildcards.sample}.sh "