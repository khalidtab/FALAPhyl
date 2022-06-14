rule alpha_div_calc: # Provides per sample alpha calculation
   version: "1.0"
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      alphadiv=report("data/alpha_div/calc_{sample}–{alpha}.txt",
      caption="../report/alpha_calc.rst",
      category="Alpha diversity",
      subcategory="{alpha}",
      labels={
              "Data type": "Calculations - text file",
              "Method": "{alpha}"
              })
   message: "Alpha diversity - {wildcards.alpha}: Calculating alpha diversity for {wildcards.sample}"
   shell:
      "Rscript --vanilla ./workflow/scripts/alpha_generate.R -i {input} -o {output.alphadiv} -a {wildcards.alpha} "


rule alpha_div_plot: # Calculates alpha diversity in each group, and outputs a PDF plot, and calculates the alpha diversity statistical analysis using nonparametric methods
   version: "1.0"
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      map="data/{sample}.txt",
      alphadiv="data/alpha_div/calc_{sample}–{alpha}.txt"
   output:
      svg=report("data/plots/alpha_div_{sample}/{group}–{alpha}.svg",
      caption="../report/alpha_plot.rst",
      category="Alpha diversity",
      subcategory="{alpha}",
      labels={
              "Data type": "Violin plot",
              "Method": "{alpha}",
              "Grouping category": "{group}"
              })
   params: 
      color=config["color"][0],
      width=config["width"][0],
      height=config["height"][0]
   message: "Alpha diversity - {wildcards.alpha}: Plotting variable {wildcards.group} for {wildcards.sample}"
   shell:
      "mkdir -p tmp data/plots/alpha_div_{wildcards.sample} &&"
      "Rscript --vanilla ./workflow/scripts/alpha_plot.R -i {input.alphadiv} -m {input.map} -c {wildcards.group}{params.color} -g {wildcards.group} -o {output.svg} -x {params.width} -y {params.height} "


rule alpha_div_stats: # Provides per sample alpha calculation
   version: "1.0"
   conda:
      "../../workflow/envs/ggrepel.yaml"
   input:
      "data/alpha_div/calc_{sample}–{alpha}.txt"
   output:
      alphadiv=report("data/alpha_div/stats_{sample}–{group}–{alpha}.txt",
      caption="../report/alpha_stats.rst",
      category="Alpha diversity",
      subcategory="{alpha}",
      labels={
              "Data type": "Nonparametric statistical testing - text file",
              "Method": "{alpha}",
              "Grouping category": "{group}"
              })
   message: "Alpha diversity - {wildcards.alpha}: Performing non-parametric testing on  {wildcards.sample}'s variable {wildcards.group}"
   shell:
      "Rscript --vanilla ./workflow/scripts/alpha_div_nonparam.R -i data/alpha_div/calc_{wildcards.sample}–{wildcards.alpha}.txt -o {output.alphadiv} -m data/{wildcards.sample}.txt -g {wildcards.group} 2>> /dev/null"


rule alpha_div: # Final step in alpha diversity calculations, to clean up temporary files
   version: "1.0"
   input:
      svg="data/plots/alpha_div_{sample}/{group}–{alpha}.svg",
      stats="data/alpha_div/stats_{sample}–{group}–{alpha}.txt"
   output:
      touch("tmp/expand–{sample}–{alpha}–{group}.txt")
   message: "Alpha diversity: Cleaning up…"
   shell:
      "find data/plots/alpha_* -empty -type d -delete &&"
      "find data/alpha_div/ -empty -type d -delete "


rule alpha:
  input:
   expand("tmp/expand–{sample}–{alpha}–{group}.txt", sample=config["mysample"], alpha=config["alpha"], group=config["group"])