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
      "; done' > temp/ANOSIM_{wildcards.sample}.sh &&"
      "chmod +x temp/ANOSIM_{wildcards.sample}.sh &&"
      "bash temp/ANOSIM_{wildcards.sample}.sh"

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
      "; done' > temp/PERMDISP_{wildcards.sample}.sh &&"
      "chmod +x temp/PERMDISP_{wildcards.sample}.sh &&"
      "bash temp/PERMDISP_{wildcards.sample}.sh "

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
      "; done ; done' > temp/ADONIS_{wildcards.sample}.sh &&"
      "chmod +x temp/ADONIS_{wildcards.sample}.sh &&"
      "bash temp/ADONIS_{wildcards.sample}.sh"


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
      "mkdir -p data/logs &&"
      "echo 'for x in {params.dist}; do weasyprint data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x/index.html data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x.pdf 2>> data/logs/PDF_{wildcards.sample}_ANOSIM.log && rm -rf data/distance/ANOSIM/ANOSIM_{wildcards.sample}_$x"
      "; done' > temp/PDF_ANOSIM_{wildcards.sample}.sh &&"
      "chmod +x temp/PDF_ANOSIM_{wildcards.sample}.sh &&"
      "bash temp/PDF_ANOSIM_{wildcards.sample}.sh "


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
      "mkdir -p data/logs &&"
      "echo 'for x in {params.dist}; do weasyprint data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x/index.html data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x.pdf 2>> data/logs/PDF_{wildcards.sample}_PERMDISP.log && rm -rf data/distance/PERMDISP/PERMDISP_{wildcards.sample}_$x"
      "; done' > temp/PDF_PERMDISP_{wildcards.sample}.sh &&"
      "chmod +x temp/PDF_PERMDISP_{wildcards.sample}.sh &&"
      "bash temp/PDF_PERMDISP_{wildcards.sample}.sh "


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
      "mkdir -p data/logs &&"
      "echo 'for x in {params.dist}; do for y in {params.formula}; do weasyprint data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x/index.html data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x.pdf 2>> data/logs/PDF_{wildcards.sample}_ADONIS.log && rm -rf data/distance/ADONIS/ADONIS_{wildcards.sample}_$y+$x"
      "; done ; done' > temp/PDF_ADONIS_{wildcards.sample}.sh &&"
      "chmod +x temp/PDF_ADONIS_{wildcards.sample}.sh &&"
      "bash temp/PDF_ADONIS_{wildcards.sample}.sh "
