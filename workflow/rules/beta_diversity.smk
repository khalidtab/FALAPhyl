def get_mem_mb(wildcards, attempt):
    return attempt * 1000

rule beta_div: # Calculate distances between samples based on the chosen beta diversity choices in the input file
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   resources:
      mem_mb=get_mem_mb
   params: 
      dist="{dist}"
   log:
      "data/logs/betaDiv–{sample}–{dist}.txt"
   output:
      tsv=report("data/distance/beta_div/{sample}–{dist}.tsv",
      caption="../report/beta_matrix.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Data type": "Distance Matrix",
              "Distance type": "{dist}"})
   message: "Beta diversity: Calculating beta dissimilarity {wildcards.dist} for {wildcards.sample}"
   shell: 
      "Rscript --vanilla ./workflow/scripts/beta_diversity.R -i {input} -o {output} -d {wildcards.dist} > {log}"

rule nmds_: # Plots the beta diversity distances using the Non-Metric Dimensional Scaling (NMDS) algorithm
   conda: "../../workflow/envs/ggrepel.yaml"
   input:
      rules.beta_div.output.tsv
   output: 
      mysvg=report("data/plots/betaDiv_{sample}/NMDS–{dist}–{group}.svg",
      caption="../report/beta_nmds.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "NMDS",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      mysvgpdf=report("data/plots/betaDiv_{sample}/NMDS–{dist}–{group}.svgnoNameswithProbDF.svg",
      caption="../report/beta_nmds.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "NMDS with probability density functions",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      mysvgnames=report("data/plots/betaDiv_{sample}/NMDS–{dist}–{group}.svgwithnames.svg",
      caption="../report/beta_nmds.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "NMDS with sample names",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      mysvgnamespdf= report("data/plots/betaDiv_{sample}/NMDS–{dist}–{group}.svgwithnamesProbDF.svg",
      caption="../report/beta_nmds.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "NMDS with sample names and probability density functions",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"})
   params:
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"]),
      color=config["color"][0],
      width=config["width"][0],
      height=config["height"][0]
   log:
      "data/logs/NMDS_{sample}–{dist}–{group}.log"
   message: "Beta diversity - {wildcards.dist}: Plotting NMDS for variable {wildcards.group} in {wildcards.sample}"
   shell:
      "Rscript --vanilla ./workflow/scripts/NMDS.R -i {input} -o {output.mysvg} -m data/{wildcards.sample}.txt -g {wildcards.group} -c {wildcards.group}{params.color} -x {params.width} -y {params.height} > {log} 2>&1 "

rule pcoa_: # Plots the beta diversity distances using the Principal Coordinates Analysis (PCoA) algorithm
   conda:
      "../../workflow/envs/ggrepel.yaml"
   input:
      betaDiv=rules.beta_div.output
   params: 
      dist="{dist}",
      group="{group}",
      color=config["color"][0],
      width=config["width"][0],
      height=config["height"][0]
   log:
      "data/logs/PCoA–{sample}–{dist}–{group}.txt"
   output:
      mysvg=report("data/plots/betaDiv_{sample}/PCoA–{dist}–{group}.svg",
      caption="../report/beta_pcoa.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "PCoA",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      mysvgpdf=report("data/plots/betaDiv_{sample}/PCoA–{dist}–{group}.svgnoNameswProbDF.svg",
      caption="../report/beta_pcoa.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "PCoA with probability density functions",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      mysvgnames=report("data/plots/betaDiv_{sample}/PCoA–{dist}–{group}.svgwithnames.svg",
      caption="../report/beta_pcoa.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "PCoA Plot with sample names",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      mysvgnamespdf=report("data/plots/betaDiv_{sample}/PCoA–{dist}–{group}.svgwithnamesprobDF.svg",
      caption="../report/beta_pcoa.rst", 
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "PCoA Plot with sample names and probability density functions",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"})
   message: "Beta diversity - {wildcards.dist}: Plotting PCoA for variable {wildcards.group} in {wildcards.sample}"
   shell: 
      "Rscript --vanilla ./workflow/scripts/PCoA.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -c {wildcards.group}{params.color} -o {output.mysvg} -x {params.width} -y {params.height} > {log} 2>&1 "


rule adonis_: # Calculates whether the two groups have similar dispersions (variances) to their centroid
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      rules.beta_div.output
   output:
      myresults=report("data/distance/ADONIS/{sample}/adonis–{dist}–{group}.txt",
      caption="../report/beta_adonis.rst",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "ADONIS",
              "Data type": "Text file",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      mypairwise=report("data/distance/ADONIS/{sample}/adonis–{dist}–{group}.txt_pairwise.txt",
      caption="../report/beta_adonis.rst",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "ADONIS - Pairwise with FDR correction",
              "Data type": "Text file",
              "Distance type": "{dist}",
              "Grouping category": "{group}"})
   resources:
      mem_mb=get_mem_mb
   threads: workflow.cores * 0.5
   log:
      "data/logs/ADONIS–{sample}–{dist}–{group}.txt"
   params: 
      dist="{dist}",
      group="{group}",
      color=config["color"][0],
      test="adonis",
      width=config["width"][0],
      height=config["height"][0]
   message: "Beta diversity - {wildcards.dist}: Calculating ADONIS for variable {wildcards.group} in {wildcards.sample}"
   shell:
      " Rscript --vanilla ./workflow/scripts/adonis_anosim_betadisper.R -i {input} -o {output.myresults} -m data/{wildcards.sample}.txt -g {wildcards.group} -c {params.color} -t {params.test} -x {params.width} -y {params.height} > {log} 2>&1 "

use rule adonis_ as anosim_ with:
   output:
      myresults=report("data/distance/ANOSIM/{sample}/anosim–{dist}–{group}.txt",
      caption="../report/beta_anosim.rst",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "ANOSIM",
              "Data type": "Text file",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      mypairwise=report("data/distance/ANOSIM/{sample}/anosim–{dist}–{group}_pairwise.txt",
      caption="../report/beta_anosim.rst",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "ANOSIM - Pairwise with FDR correction",
              "Data type": "Text file",
              "Distance type": "{dist}",
              "Grouping category": "{group}"})
   log: 
      "data/logs/ANOSIM–{sample}–{dist}–{group}.txt"
   params: 
      dist="{dist}",
      group="{group}",
      color=config["color"][0],
      test="anosim",
      width=config["width"][0],
      height=config["height"][0]
   message: "Beta diversity - {wildcards.dist}: Calculating ANOSIM for variable {wildcards.group} in {wildcards.sample}"


rule permdisp_: # Calculates whether the two groups have similar dispersions (variances) to their centroid
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      rules.beta_div.output
   output:
      myresults=report("data/distance/PERMDISP/{sample}/betadisper–{dist}–{group}.txt",
      caption="../report/beta_permdisp.rst",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "Beta dispersion and pairwise comparisons with FDR correction",
              "Data type": "Text file",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      myPCoA=report("data/plots/PCoA_betadispersion–{sample}–{dist}–{group}.svg",
      caption="../report/beta_permdisp.rst",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "PCoA Plot with standard deviation ellipsis of the dispersion from group centroid",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"}),
      myBoxplot=report("data/plots/boxplot_betadispersion–{sample}–{dist}–{group}.svg",
      caption="../report/beta_permdisp.rst",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "Boxplot of the beta dispersion",
              "Data type": "Plot - SVG",
              "Distance type": "{dist}",
              "Grouping category": "{group}"})
   resources:
      mem_mb=get_mem_mb
   threads: workflow.cores * 0.5
   log: 
      "data/logs/PERMDISP–{sample}–{dist}–{group}.txt"
   params: 
      dist="{dist}",
      group="{group}",
      color=config["color"][0],
      test="betadisper",
      width=config["width"][0],
      height=config["height"][0]
   message: "Beta diversity - {wildcards.dist}: Calculating beta dispersion for variable {wildcards.group} in {wildcards.sample}"
   shell:
      " Rscript --vanilla ./workflow/scripts/adonis_anosim_betadisper.R -i {input} -o {output.myresults} -m data/{wildcards.sample}.txt -p {output.myPCoA} -b {output.myBoxplot} -g {wildcards.group} -c {params.color} -t {params.test} -x {params.width} -y {params.height} > {log} 2>&1 "

rule dunns_: # Plots the beta diversity distances using the Non-Metric Dimensional Scaling (NMDS) algorithm
   conda: "../../workflow/envs/dunn.yaml"
   input:
      rules.beta_div.output.tsv
   output: 
      report("data/distance/PerGroupDistCompare/{sample}/Dunns–{dist}–{group}.txt",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "Dunn's multiple comparisons of the distances between the groups",
              "Data type": "Text file",
              "Distance type": "{dist}",
              "Grouping category": "{group}"})
   log:
      "data/logs/Dunn_{sample}–{dist}–{group}.log"
   message: "Beta diversity - {wildcards.dist}: Calculating distances between groups based on Dunn's multiple comaparison for variable {wildcards.group} in {wildcards.sample}"
   shell:
      "Rscript --vanilla ./workflow/scripts/Dunn.R -i {input} -o {output} -m data/{wildcards.sample}.txt -c {wildcards.group} > {log} 2>&1 "

rule betadisper:
   input:
      expand("data/distance/PERMDISP/{sample}/betadisper–{dist}–{group}.txt", sample=config["mysample"], dist=config["distances"], group=config["group"])

rule anosim:
   input:
      expand("data/distance/ANOSIM/{sample}/anosim–{dist}–{group}.txt",       sample=config["mysample"], dist=config["distances"], group=config["group"])

rule adonis:
   input:
      expand("data/distance/ADONIS/{sample}/adonis–{dist}–{group}.txt",       sample=config["mysample"], dist=config["distances"], group=config["group"])

rule nmds:
   input:
      expand("data/plots/betaDiv_{sample}/NMDS–{dist}–{group}.svg",           sample=config["mysample"], dist=config["distances"], group=config["group"])

rule pcoa:
   input:
      expand("data/plots/betaDiv_{sample}/PCoA–{dist}–{group}.svg",           sample=config["mysample"], dist=config["distances"], group=config["group"])

rule dunn:
   input:
      expand("data/distance/PerGroupDistCompare/{sample}/Dunns–{dist}–{group}.txt",           sample=config["mysample"], dist=config["distances"], group=config["group"])

rule beta:
   input:
      expand("data/distance/PERMDISP/{sample}/betadisper–{dist}–{group}.txt", sample=config["mysample"], dist=config["distances"], group=config["group"]),
      expand("data/distance/ANOSIM/{sample}/anosim–{dist}–{group}.txt",       sample=config["mysample"], dist=config["distances"], group=config["group"]),
      expand("data/distance/ADONIS/{sample}/adonis–{dist}–{group}.txt",       sample=config["mysample"], dist=config["distances"], group=config["group"]),
      expand("data/plots/betaDiv_{sample}/NMDS–{dist}–{group}.svg",           sample=config["mysample"], dist=config["distances"], group=config["group"]),
      expand("data/plots/betaDiv_{sample}/PCoA–{dist}–{group}.svg",           sample=config["mysample"], dist=config["distances"], group=config["group"]),
      expand("data/distance/PerGroupDistCompare/{sample}/Dunns–{dist}–{group}.txt",           sample=config["mysample"], dist=config["distances"], group=config["group"])
      
