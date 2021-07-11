rule make_tree: #If a tree is not provided, a tree is generated from the feature table using the Ward hierarchical clustering which is necessary for using PhILR.
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      ancient("data/biom/{sample}.biom")
   output:
      "data/tree/{sample}.tre"
   message: "Tree for {wildcards.sample} does not exist. Generating one based on Ward hierarchical clustering."
   shell:
      "   mkdir -p data/tsv &&"
      "   biom convert -i {input} -o data/tsv/{wildcards.sample}_for_tree.tsv --to-tsv &&"
      '   biom convert -i data/tsv/{wildcards.sample}_for_tree.tsv -o data/biom/{wildcards.sample}_for_tree.biom --to-json --table-type="OTU table" &&'
      "   qiime tools import --input-path data/biom/{wildcards.sample}_for_tree.biom --output-path data/biom/{wildcards.sample}_for_tree.qza --type FeatureTable[Frequency] --input-format BIOMV100Format >/dev/null &&"
      "   qiime gneiss correlation-clustering --i-table data/biom/{wildcards.sample}_for_tree.qza --o-clustering data/tree/{wildcards.sample}.qza &&"
      "   qiime tools export --input-path data/tree/{wildcards.sample}.qza --output-path data/tree/tree_{wildcards.sample}  >/dev/null &&"
      "   mv data/tree/tree_{wildcards.sample}/tree.nwk {output} &&"
      "   rm data/tree/{wildcards.sample}.qza data/biom/{wildcards.sample}_for_tree.qza data/tsv/{wildcards.sample}_for_tree.tsv &&"
      "   rm -rf data/tree/tree_{wildcards.sample}"



rule philr_distance: # Calculates PhILR distances between different samples, without taking into account phylogenic distances (this can be changed by modifying the philr.R script)
   version: "1.0"
   conda:
      "../../workflow/envs/philr.yaml"
   input:
      biom = "data/biom/{sample}.biom",
      tree = "data/tree/{sample}.tre",
      map = "data/map/{sample}.txt"
   output: "data/distance/beta_div/{sample}+philr.tsv"
   message: "Creating PhILR distances for {wildcards.sample}"
   log: "data/logs/philr_{sample}.log"
   resources:
      mem_mb=10000
   shell:
      "mkdir -p data/distance/beta_div/ data/logs &&"
      "Rscript --vanilla ./workflow/scripts/philr.R -i {input.biom} -t {input.tree} -m {input.map} -o {output} 2>> data/logs/philr_{wildcards.sample}.log"
      

rule philr_processing_in_qiime2: # Importing PhILR distances to qiime2 for calculating PCoA coordinates.
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      "data/distance/beta_div/{sample}+philr.tsv"
   output: 
      philr_pcoa_tsv="data/distance/PCoA/PCoA_{sample}+philr.tsv",
      philr_dist_qza=temporary("data/distance/beta_div/{sample}+philr.qza"),
      philr_pcoa_qza=temporary("data/distance/PCoA/PCoA_{sample}+philr.qza")
   message: "Importing newly created PhILR distance for analysing {wildcards.sample}, and generating distances/PCoA coordinates."
   shell:
      "mkdir -p data/distance/beta_div/ data/distance/PCoA/ &&"
      "qiime tools import --input-path {input} --output-path data/distance/beta_div/{wildcards.sample}+philr.qza --type DistanceMatrix && "
      "qiime diversity pcoa --i-distance-matrix data/distance/beta_div/{wildcards.sample}+philr.qza --o-pcoa data/distance/PCoA/PCoA_{wildcards.sample}+philr.qza  >/dev/null &&"
      "qiime tools export --input-path data/distance/PCoA/PCoA_{wildcards.sample}+philr.qza --output-path data/distance/PCoA/PCoA_{wildcards.sample}+philr  >/dev/null &&"
      "cp data/distance/PCoA/PCoA_{wildcards.sample}+philr/ordination.txt data/distance/PCoA/PCoA_{wildcards.sample}+philr.tsv && "
      "rm -rf data/distance/PCoA/PCoA_{wildcards.sample}+philr/"

rule philr_pcoa: # Plots the PhILR PCoA coordinates using the Principal Coordinates Analysis (PCoA) algorithm
   version: "1.0"
   conda:
      "../../workflow/envs/phylotoast.yaml"
   input:
      betaDiv=rules.philr_processing_in_qiime2.output.philr_pcoa_tsv,
      color=rules.get_colors.output
   params: 
      group=expand("{group}",group=config["group"])
   output:
      report(expand("data/plots/PCoA_{{sample}}+philr+{group}.svg",group=config["group"]))
   message: "Generating PhILR-PCoA plots for {wildcards.sample}"
   shell: 
      "mkdir -p data/plots &&"
      "echo 'for y in {params.group};" 
      "do xvfb-run --auto-servernum python2 workflow/scripts/PCoA.py -i data/distance/PCoA/PCoA_{wildcards.sample}+philr.tsv -m data/map/{wildcards.sample}.txt -b $y -d 2 -c data/map/color_{wildcards.sample}+$y.txt -o data/plots/PCoA_{wildcards.sample}+philr+$y.svg;  done'"
      "> tmp/SVG_PCoA_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/SVG_PCoA_philr_{wildcards.sample}.sh &&"
      "bash tmp/SVG_PCoA_philr_{wildcards.sample}.sh "
      
      
rule philr_nmds: # Plots the PhILR distances using the Non-Metric Dimensional Scaling (NMDS) algorithm
   version: "1.0"
   conda: "../../workflow/envs/NMDS_plotly.yaml"
   input:
      rules.philr_distance.output
   output: 
      json=temporary(expand("data/plots/NMDS_{{sample}}+philr+{group}.json",group=config["group"])),
      svg=report(expand("data/plots/NMDS_{{sample}}+philr+{group}_1.svg",group=config["group"]))
   params:
      group=expand("{group}",group=config["group"]),
      color=config["color"][0]
   log:
      expand("data/logs/NMDS_{{sample}}+philr+{group}.log",group=config["group"])
   message: "Generating PhILR-NMDS plots for {wildcards.sample}"
   shell:
      "echo 'for w in {params.group}; do "
      "y=$(printf \"%s{params.color}\" $w) && "
      "Rscript --vanilla ./workflow/scripts/NMDS.R -i data/distance/beta_div/{wildcards.sample}+philr.tsv -o data/plots/NMDS_{wildcards.sample}+philr+$w.json -m data/map/{wildcards.sample}.txt -g $w -c $y > data/logs/NMDS_{wildcards.sample}+philr+$w.log 2>>/dev/null && "
      "xvfb-run --auto-servernum orca graph data/plots/NMDS_{wildcards.sample}+philr+$w.json -o data/plots/NMDS_{wildcards.sample}+philr+$w.svg -f svg 2>> data/logs/NMDS_{wildcards.sample}+philr+$w.log || true; done' > tmp/beta_div_NMDS_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_NMDS_philr_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_NMDS_philr_{wildcards.sample}.sh"

rule philr_nmds_hull: # Plots the PhILR distances using the Non-Metric Dimensional Scaling (NMDS) algorithm, and adds a hull around the samples
   version: "1.0"
   conda: "../../workflow/envs/NMDS_plotly.yaml"
   input:
      rules.philr_distance.output
   output: 
      json=temporary(expand("data/plots/NMDS_HULL_{{sample}}+philr+{group}.json",group=config["group"])),
      svg=report(expand("data/plots/NMDS_HULL_{{sample}}+philr+{group}_1.svg",group=config["group"]))
   params:
      group=expand("{group}",group=config["group"]),
      color=config["color"][0]
   log:
      expand("data/logs/NMDS_HULL_{{sample}}+philr+{group}.log",group=config["group"])
   message: "Generating PhILR NMDS hull plots for {wildcards.sample}"
   shell:
      "echo 'for w in {params.group}; do "
      "y=$(printf \"%s{params.color}\" $w) && "
      " Rscript --vanilla ./workflow/scripts/NMDS_hull.R -i data/distance/beta_div/{wildcards.sample}+philr.tsv -o data/plots/NMDS_HULL_{wildcards.sample}+philr+$w.json -m data/map/{wildcards.sample}.txt -g $w -c $y > data/logs/NMDS_HULL_{wildcards.sample}+philr+$w.log 2>>/dev/null && "
      "xvfb-run --auto-servernum orca graph data/plots/NMDS_HULL_{wildcards.sample}+philr+$w.json -o data/plots/NMDS_HULL_{wildcards.sample}+philr+$w.svg -f svg 2>> data/logs/NMDS_HULL_{wildcards.sample}+philr+$w.log || true; done' > tmp/beta_div_NMDS_HULL_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/beta_div_NMDS_HULL_philr_{wildcards.sample}.sh &&"
      "bash tmp/beta_div_NMDS_HULL_philr_{wildcards.sample}.sh"




rule philr_anosim: # Calculates whether the intra-group PhILR variances is sig different from intergroup variances
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.philr_processing_in_qiime2.output.philr_dist_qza
   output:
            temporary(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+philr.qzv", group=config["group"]))
   params:
      group=expand("{group}",group=config["group"])
   message: "Calculating PhILR ANOSIM for {wildcards.sample}"
   shell:
      "mkdir -p data/distance/ANOSIM &&"
      "echo 'for y in {params.group}; do qiime diversity beta-group-significance --i-distance-matrix {input} --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column $y --p-method anosim --p-pairwise --o-visualization data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+philr.qzv >/dev/null && "
      "qiime tools export --input-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+philr.qzv --output-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+philr >/dev/null"
      "; done' > tmp/ANOSIM_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/ANOSIM_philr_{wildcards.sample}.sh &&"
      "bash tmp/ANOSIM_philr_{wildcards.sample}.sh"

rule make_philr_anosim_PDFs: # Creates PDFs out of the HTML file for easier handling
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.philr_anosim.output
   output: 
      anosim=report(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+philr.pdf", group=config["group"]))
   log: "data/logs/PDF_{sample}+philr+ANOSIM.log"
   params: 
      group=expand("{group}",group=config["group"])
   message: "Creating PhILR ANOSIM PDFs for {wildcards.sample}"
   shell:
      "mkdir -p data/logs &&"
      "echo 'for y in {params.group}; do weasyprint data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+philr/index.html data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+philr.pdf 2>> data/logs/PDF_{wildcards.sample}+philr+ANOSIM.log && rm -rf data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+philr"
      "; done ' > tmp/PDF_ANOSIM_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/PDF_ANOSIM_philr_{wildcards.sample}.sh &&"
      "bash tmp/PDF_ANOSIM_philr_{wildcards.sample}.sh "



rule philr_permdisp: # Calculates whether the two groups have similar dispersions (variances) of their PhILR distances to their centroid
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.philr_processing_in_qiime2.output.philr_dist_qza
   output:
      temporary(expand("data/distance/PERMDISP/PERMDISP_{{sample}}+{group}+philr.qzv", group=config["group"]))
   params: 
      group=expand("{group}",group=config["group"])
   message: "Calculating PhILR PERMDISP for {wildcards.sample}"
   shell:
      "mkdir -p data/distance/PERMDISP &&"
      "echo 'for y in {params.group}; do qiime diversity beta-group-significance --i-distance-matrix {input} --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column $y --p-method permdisp --p-pairwise --o-visualization data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+philr.qzv >/dev/null && "
      "qiime tools export --input-path data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+philr.qzv --output-path data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+philr >/dev/null"
      "; done' > tmp/PERMDISP_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/PERMDISP_philr_{wildcards.sample}.sh &&"
      "bash tmp/PERMDISP_philr_{wildcards.sample}.sh"

rule make_philr_permdisp_PDFs: # Creates PDFs out of the HTML file for easier handling
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.philr_permdisp.output
   output: 
      permdisp=report(expand("data/distance/PERMDISP/PERMDISP_{{sample}}+{group}+philr.pdf", group=config["group"]))
   log: "data/logs/PDF_{sample}+philr+PERMDISP.log"
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"])
   message: "Creating PhILR PERMDISP PDFs for {wildcards.sample}"
   shell:
      "mkdir -p data/logs &&"
      "echo 'for y in {params.group}; do weasyprint data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+philr/index.html data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+philr.pdf 2>> data/logs/PDF_{wildcards.sample}+philr+PERMDISP.log && rm -rf data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+philr"
      "; done' > tmp/PDF_PERMDISP_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/PDF_PERMDISP_philr_{wildcards.sample}.sh &&"
      "bash tmp/PDF_PERMDISP_philr_{wildcards.sample}.sh "

rule philr_adonis: # Calculates whether the two groups have different centroids, susceptible to the groups having different dispersions. Therefore, interpret along with PERDISP (also known as betadisper)
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.philr_processing_in_qiime2.output.philr_dist_qza
   output:
            temporary(expand("data/distance/ADONIS/ADONIS_{{sample}}+{group}+philr.qzv", group=config["group"]))
   params: 
      group=expand("{group}",group=config["group"])
   message: "Calculating PhILR ADONIS for {wildcards.sample}"
   shell:
      "mkdir -p data/distance/ADONIS &&"
      "echo 'for y in {params.group}; do qiime diversity adonis --i-distance-matrix {input} --m-metadata-file data/map/{wildcards.sample}.txt --p-formula \"$y\" --o-visualization data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+philr.qzv >/dev/null && "
      "qiime tools export --input-path data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+philr.qzv --output-path data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+philr >/dev/null"
      "; done' > tmp/ADONIS_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/ADONIS_philr_{wildcards.sample}.sh &&"
      "bash tmp/ADONIS_philr_{wildcards.sample}.sh"


rule make_philr_adonis_PDFs: # Creates PDFs out of the HTML file for easier handling
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.philr_adonis.output
   output: 
      adonis=report(expand("data/distance/ADONIS/ADONIS_{{sample}}+{group}+philr.pdf", group=config["group"]))
   log: "data/logs/PDF_{sample}+ADONIS.log"
   params: 
      group=expand("{group}",group=config["group"])
   message: "Creating PhILR ADONIS PDFs for {wildcards.sample}"
   shell:
      "mkdir -p data/logs &&"
      "echo 'for y in {params.group}; do weasyprint data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+philr/index.html data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+philr.pdf 2>> data/logs/PDF_{wildcards.sample}+philr+ADONIS.log && rm -rf data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+philr"
      "; done' > tmp/PDF_ADONIS_philr_{wildcards.sample}.sh &&"
      "chmod +x tmp/PDF_ADONIS_philr_{wildcards.sample}.sh &&"
      "bash tmp/PDF_ADONIS_philr_{wildcards.sample}.sh "

rule compositional_betadiv: # final step of compositional analysis
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      adonis=rules.make_philr_adonis_PDFs.output,
      permdisp=rules.make_philr_permdisp_PDFs.output,
      anosim=rules.make_philr_anosim_PDFs.output,
      nmds=rules.philr_nmds.output.svg,
      nmdshull=rules.philr_nmds_hull.output.svg,
      pcoa=rules.philr_pcoa.output
   output: 
      touch(temporary("tmp/compositional_{sample}.final"))
   message: "Finalizing compositional analysis for {wildcards.sample}"
   shell:
      "echo 'Compositional analysis of {wildcards.sample} is done.'"
