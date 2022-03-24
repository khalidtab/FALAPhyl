rule alpha_div_calc: # Provides per sample alpha calculation
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      qza=rules.biom_to_qza.output.qza
   output:
      alphadiv=report(expand("data/alpha_div/calc_{{sample}}+{alpha}.txt",alpha=config["alpha"]))
   params: 
      alpha=expand("{alpha}",alpha=config["alpha"])
   message: "Calculating alpha diversity for {wildcards.sample}"
   shell:
      "mkdir -p tmp &&"
      "echo 'for x in {params.alpha}; do mkdir -p data/alpha_div/ && "
      "qiime diversity alpha --i-table {input.qza} --p-metric $x --o-alpha-diversity data/alpha_div/{wildcards.sample}+$x.qza && "
      "qiime tools export --input-path data/alpha_div/{wildcards.sample}+$x.qza --output-path data/alpha_div/calc_{wildcards.sample}+$x  >/dev/null && "
      "rm data/alpha_div/{wildcards.sample}+$x.qza && "
      "mv data/alpha_div/calc_{wildcards.sample}+$x/alpha-diversity.tsv data/alpha_div/calc_{wildcards.sample}+$x.txt "
      "; done' > tmp/alpha_div_{wildcards.sample}.sh && "
      "chmod +x tmp/alpha_div_{wildcards.sample}.sh && "
      "bash tmp/alpha_div_{wildcards.sample}.sh"

rule alpha_div_plot: # Calculates alpha diversity in each group, and outputs a PDF plot, and calculates the alpha diversity statistical analysis using nonparametric methods
   version: "1.0"
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      map="data/map/{sample}.txt",
      alphadiv=rules.alpha_div_calc.output
   output:
      svg=report(expand("data/plots/alpha_div_{{sample}}+{group}+{alpha}.svg",alpha=config["alpha"],group=config["group"]))
   params: 
      alpha=expand("{alpha}",alpha=config["alpha"]),
      group=expand("{group}",group=config["group"]),
      color=config["color"][0]
   message: "Plotting alpha diversity on {wildcards.sample}"
   shell:
      "mkdir -p tmp &&"
      "echo 'for x in {params.alpha}; do for y in {params.group}; do "
      " w=$(printf \"%s{params.color}\" $y) && "
      "mkdir -p data/plots/alpha_div_{wildcards.sample}+$y+$x && "
      "Rscript --vanilla ./workflow/scripts/alphaDiv.R -i data/alpha_div/calc_{wildcards.sample}+$x.txt -m {input.map} -c $w -g $y -o data/plots/alpha_div_{wildcards.sample}+$y+$x.svg "
      "; done ; done' > tmp/alpha_div_plot_{wildcards.sample}.sh && "
      "chmod +x tmp/alpha_div_plot_{wildcards.sample}.sh && "
      "bash tmp/alpha_div_plot_{wildcards.sample}.sh"



rule alpha_div_stats: # Provides per sample alpha calculation
   version: "1.0"
   conda:
      "../../workflow/envs/ggrepel.yaml"
   input:
      qza=rules.alpha_div_calc.output
   output:
      alphadiv=report(expand("data/alpha_div/stats_{{sample}}+{group}+{alpha}.txt",alpha=config["alpha"],group=config["group"]))
   params: 
      alpha=expand("{alpha}",alpha=config["alpha"]),
      group=expand("{group}",group=config["group"])
   message: "Performing non-parametric testing on alpha diversity of {wildcards.sample}"
   shell:
      "echo 'for x in {params.alpha}; do for y in {params.group}; do mkdir -p data/alpha_div/ && "
      "Rscript --vanilla ./workflow/scripts/alpha_div_nonparam.R -i data/alpha_div/calc_{wildcards.sample}+$x.txt -o data/alpha_div/stats_{wildcards.sample}+$y+$x.txt -m data/map/{wildcards.sample}.txt -g $y 2>> /dev/null"
      "; done; done' > tmp/stats_alpha_div_{wildcards.sample}.sh && "
      "chmod +x tmp/stats_alpha_div_{wildcards.sample}.sh && "
      "bash tmp/stats_alpha_div_{wildcards.sample}.sh"





rule alpha_div: # Final step in alpha diversity calculations, to clean up temporary files
   version: "1.0"
   input:
      svg=rules.alpha_div_plot.output.svg,
      stats=rules.alpha_div_stats.output.alphadiv
   output:
      temporary(touch("tmp/{sample}_alpha_div_all_done.tmp"))
   message: "Cleaning up from alpha diversity"
   shell:
      "find data/plots/alpha_* -empty -type d -delete &&"
      "find data/alpha_div/ -empty -type d -delete "