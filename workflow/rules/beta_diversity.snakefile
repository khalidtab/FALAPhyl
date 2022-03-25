rule beta_div: # Calculate distances between samples based on the chosen beta diversity choices in the input file
   version: "2.0"
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   params: 
      dist=expand("{dist}",dist=config["distances"])
   output:
      tsv=expand("data/distance/beta_div/{{sample}}+{dist}.tsv",dist=config["distances"])
   message: "Calculating beta diversity for {wildcards.sample}"
   shell: 
      "mkdir -p data/distance/beta_div &&"
      "echo 'for x in {params.dist}; do "
      "Rscript --vanilla ./workflow/scripts/beta_diversity.R -i {input} -o /data/distance/beta_div/{wildcards.sample}+$x.tsv -d $x"
      "; done' > tmp/beta_div_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_{wildcards.sample}.sh"

rule nmds: # Plots the beta diversity distances using the Non-Metric Dimensional Scaling (NMDS) algorithm
   version: "1.0"
   conda: "../../workflow/envs/ggrepel.yaml"
   input:
      rules.beta_div.output.tsv
   output: 
      report(expand("data/plots/NMDS_{{sample}}+{dist}+{group}.svg",dist=config["distances"],group=config["group"]))
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
      "Rscript --vanilla ./workflow/scripts/NMDS.R -i data/distance/beta_div/{wildcards.sample}+$x.tsv -o data/plots/NMDS_{wildcards.sample}+$x+$w.svg -m data/map/{wildcards.sample}.txt -g $w -c $y > data/logs/NMDS_{wildcards.sample}+$x+$w.log 2>>/dev/null; done; done' > tmp/SVG_NMDS_{wildcards.sample}.sh &&"
      "chmod +x tmp/SVG_NMDS_{wildcards.sample}.sh &&"
      "bash tmp/SVG_NMDS_{wildcards.sample}.sh"

rule pcoa: # Plots the beta diversity distances using the Principal Coordinates Analysis (PCoA) algorithm
   version: "1.0"
   conda:
      "../../workflow/envs/ggrepel.yaml"
   input:
      betaDiv=rules.beta_div.output
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"]),
      color=config["color"][0]
   output:
      report(expand("data/plots/PCoA_{{sample}}+{dist}+{group}.svg",dist=config["distances"],group=config["group"]))
   message: "Generating PCoA plots for {wildcards.sample}"
   shell: 
      "mkdir -p data/plots &&"
      "echo 'for x in {params.dist}; do for y in {params.group}; do "
      "z=$(printf \"%s{params.color}\" $y) && " 
      "Rscript --vanilla ./workflow/scripts/PCoA.R -i data/distance/PCoA/PCoA_{wildcards.sample}+$x.tsv -m data/map/{wildcards.sample}.txt -g $y -c $z -o data/plots/PCoA_{wildcards.sample}+$x+$y.svg ; done ; done'"
      "> tmp/SVG_PCoA_{wildcards.sample}.sh &&"
      "chmod +x tmp/SVG_PCoA_{wildcards.sample}.sh &&"
      "bash tmp/SVG_PCoA_{wildcards.sample}.sh "


rule adonis: # Calculates whether the two groups have similar dispersions (variances) to their centroid
   version: "1.0"
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      rules.beta_div.output
   output:
      myresults=expand("data/distance/ADONIS/{{sample}}+{dist}+{group}_adonis.txt", dist=config["distances"], group=config["group"]),
      mysh="tmp/ADONIS_{sample}.sh",
      myfolder = directory("data/distance/ADONIS/{sample}")
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"]),
      color=config["color"][0]
   message: "Calculating ADONIS for {wildcards.sample}"
   shell:
      "mkdir -p {output.myfolder} &&"
      "echo 'for x in {params.dist}; do for y in {params.group}; do"
      "Rscript --vanilla ./workflow/scripts/adonis_anosim_betadisper.R -i data/distance/beta_div/{wildcards.sample}+$x.tsv -o {output.myfolder}{wildcards.sample}+$x+$y -m data/map/{wildcards.sample}.txt -p data/plots/{wildcards.sample}+$x+$y -g $y -c {wildcards.color} -t adonis "
      "; done ; done' > {output.mysh} &&"
      "chmod +x tmp/ADONIS_{output.mysh}.sh &&"
      "bash tmp/ADONIS_{output.mysh}.sh"

use rule adonis as anosim with:
   output:
      myresults=expand("data/distance/ANOSIM/{{sample}}+{dist}+{group}_anosim.txt", dist=config["distances"], group=config["group"]),
      mysh="tmp/ANOSIM_{sample}.sh",
      myfolder = directory("data/distance/ANOSIM/{sample}")
   message: "Calculating ANOSIM for {wildcards.sample}"


use rule adonis as permdisp with: # Calculates whether the two groups have similar dispersions (variances) to their centroid
   output:
      myresults=expand("data/distance/PERMDISP/{{sample}}+{dist}+{group}_adonis.txt", dist=config["distances"], group=config["group"]),
      mysh="tmp/PERMDISP_{sample}.sh",
      myfolder = directory("data/distance/PERMDISP/{sample}")
   message: "Calculating PERMDISP for {wildcards.sample}"
