rule beta_div: #Opportunity for parallelization by writing all commands to a text file, then funnelling them to gnu parallel
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
      "mkdir -p data/distance/beta_div data/distance/PCoA data/scripts &&"
      "echo 'for x in {params.dist}; do qiime diversity beta --i-table data/biom/{wildcards.sample}.qza --o-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --p-metric $x &&"
      "qiime diversity pcoa --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --o-pcoa data/distance/PCoA/PCoA_{wildcards.sample}+$x.qza &&"
      "qiime tools export --input-path data/distance/beta_div/{wildcards.sample}+$x.qza --output-path data/distance/beta_div/{wildcards.sample}+$x &&"
      "qiime tools export --input-path data/distance/PCoA/PCoA_{wildcards.sample}+$x.qza --output-path data/distance/PCoA/PCoA_{wildcards.sample}+$x &&"
      "cp data/distance/beta_div/{wildcards.sample}+$x/distance-matrix.tsv data/distance/beta_div/{wildcards.sample}+$x.tsv &&"
      "cp data/distance/PCoA/PCoA_{wildcards.sample}+$x/ordination.txt data/distance/PCoA/PCoA_{wildcards.sample}+$x.tsv && "
      "rm -rf data/distance/PCoA/PCoA_{wildcards.sample}+$x/ && rm -rf data/distance/beta_div/{wildcards.sample}+$x/"
      "; done' > data/scripts/beta_div_PCoA_{wildcards.sample}.sh |"
      "chmod +x data/scripts/beta_div_PCoA_{wildcards.sample}.sh |"
      "bash data/scripts/beta_div_PCoA_{wildcards.sample}.sh"


rule anosim: # Calculates whether the intra-group variances is sig different from intergroup variances
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.beta_div.output.qza
   output:
      temporary(expand("data/distance/ANOSIM/ANOSIM_{{sample}}_{dist}.qzv",dist=config["distances"]))
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=config["group"][0]
   shell:
      "mkdir -p data/distance/ANOSIM &&"
      "echo 'for x in {params.dist}; do qiime diversity beta-group-significance --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column {params.group} --p-method anosim --o-visualization data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x.qzv && qiime tools export --input-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x.qzv --output-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x "
      "; done' > data/scripts/ANOSIM_{wildcards.sample}.sh |"
      "chmod +x data/scripts/ANOSIM_{wildcards.sample}.sh |"
      "bash data/scripts/ANOSIM_{wildcards.sample}.sh"

rule permdisp: # Calculates whether the two groups have similar dispersions (variances) to their centroid
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.beta_div.output.qza
   output:
      temporary(expand("data/distance/PERMDISP/PERMDISP_{{sample}}_{dist}.qzv",dist=config["distances"]))
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=config["group"][0]
   shell:
      "mkdir -p data/distance/PERMDISP &&"
      "echo 'for x in {params.dist}; do qiime diversity beta-group-significance --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column {params.group} --p-method permdisp --o-visualization data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x.qzv && qiime tools export --input-path data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x.qzv --output-path data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x"
      "; done' > data/scripts/PERMDISP_{wildcards.sample}.sh |"
      "chmod +x data/scripts/PERMDISP_{wildcards.sample}.sh |"
      "bash data/scripts/PERMDISP_{wildcards.sample}.sh "

rule adonis: # Calculates whether the two groups have different centroids, susceptible to the groups having different dispersions. Therefore, interpret along with PERDISP (also known as betadisper)
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.beta_div.output.qza
   output:
      temporary(expand("data/distance/ADONIS/ADONIS_{{sample}}_{y}+{x}.qzv", x=config["distances"], y=config["adonis"]))
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      formula=expand("{dist}",dist=config["adonis"])
   shell:
      "mkdir -p data/distance/ADONIS &&"
      "echo 'for x in {params.dist}; do for y in {params.formula}; do qiime diversity adonis --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --m-metadata-file data/map/{wildcards.sample}.txt --p-formula \"$y\" --o-visualization data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x.qzv && "
      "qiime tools export --input-path data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x.qzv --output-path data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x "
      "; done ; done' > data/scripts/ADONIS_{wildcards.sample}.sh |"
      "chmod +x data/scripts/ADONIS_{wildcards.sample}.sh |"
      "bash data/scripts/ADONIS_{wildcards.sample}.sh"


rule make_anosim_PDFs:
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.anosim.output
   output: 
      anosim=report(expand("data/distance/ANOSIM/ANOSIM_{{sample}}_{dist}.pdf",dist=config["distances"])),
   log: "data/logs/PDF_{sample}_anosim.log"
   params: 
      dist=expand("{dist}",dist=config["distances"])
   shell:
      "mkdir -p data/logs data/scripts &&"
      "echo 'for x in {params.dist}; do weasyprint data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x/index.html data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x.pdf 2>> data/logs/PDF_{wildcards.sample}_ANOSIM.log && rm -rf data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x"
      "; done' > data/scripts/PDF_ANOSIM_{wildcards.sample}.sh |"
      "chmod +x data/scripts/PDF_ANOSIM_{wildcards.sample}.sh |"
      "bash data/scripts/PDF_ANOSIM_{wildcards.sample}.sh "


rule make_permdisp_PDFs:
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.permdisp.output
   output: 
      permdisp=report(expand("data/distance/PERMDISP/PERMDISP_{{sample}}_{dist}.pdf",dist=config["distances"])),
   log: "data/logs/PDF_{sample}_permdisp.log"
   params: 
      dist=expand("{dist}",dist=config["distances"])
   shell:
      "mkdir -p data/logs data/scripts &&"
      "echo 'for x in {params.dist}; do weasyprint data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x/index.html data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x.pdf 2>> data/logs/PDF_{wildcards.sample}_PERMDISP.log && rm -rf data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x"
      "; done' > data/scripts/PDF_PERMDISP_{wildcards.sample}.sh |"
      "chmod +x data/scripts/PDF_PERMDISP_{wildcards.sample}.sh |"
      "bash data/scripts/PDF_PERMDISP_{wildcards.sample}.sh "


rule make_adonis_PDFs:
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.adonis.output
   output: 
      adonis=report(expand("data/distance/ADONIS/ADONIS_{{sample}}_{y}+{x}.pdf", x=config["distances"], y=config["adonis"]))
   log: "data/logs/PDF_{sample}_ADONIS.log"
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      formula=expand("{dist}",dist=config["adonis"])
   shell:
      "mkdir -p data/logs data/scripts &&"
      "echo 'for x in {params.dist}; do for y in {params.formula}; do weasyprint data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x/index.html data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x.pdf 2>> data/logs/PDF_{wildcards.sample}_ADONIS.log && rm -rf data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x"
      "; done ; done' > data/scripts/PDF_ADONIS_{wildcards.sample}.sh |"
      "chmod +x data/scripts/PDF_ADONIS_{wildcards.sample}.sh |"
      "bash data/scripts/PDF_ADONIS_{wildcards.sample}.sh "


rule pcoa_svg: #Opportunity for parallelization by writing all commands to a text file, then funnelling them to gnu parallel
   version: "1.0"
   conda:
      "../../workflow/envs/phylotoast.yaml"
   input:
      rules.beta_div.output.tsv_pcoa
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=config["group"][0],
      color=config["color"][0]
#   log: "data/logs/SVG_PCoA_{sample}.log"
   output:
      report(expand("data/plots/PCoA_{{sample}}+{dist}.svg",dist=config["distances"]))
   shell: 
      "mkdir -p data/plots &&"
      "echo 'for x in {params.dist}; do PCoA.py -i data/distance/PCoA/PCoA_{wildcards.sample}+$x.tsv -m data/map/{wildcards.sample}.txt -g {params.group} -d 2 -c {params.color} -o data/plots/PCoA_{wildcards.sample}+$x.svg ; done'"
      "> data/scripts/SVG_PCoA_{wildcards.sample}.sh |"
      "chmod +x data/scripts/SVG_PCoA_{wildcards.sample}.sh |"
      "bash data/scripts/SVG_PCoA_{wildcards.sample}.sh "