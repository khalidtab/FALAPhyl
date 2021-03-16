rule anosim: # Calculates whether the intra-group variances is sig different from intergroup variances
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.beta_div.output.qza
   output:
            temporary(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+{dist}.qzv", dist=config["distances"], group=config["group"]))
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"])
   message: "Calculating ANOSIM for {wildcards.sample}"
   shell:
      "mkdir -p data/distance/ANOSIM &&"
      "echo 'for x in {params.dist}; do for y in {params.group}; do qiime diversity beta-group-significance --p-pairwise --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column $y --p-method anosim --o-visualization data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+$x.qzv >/dev/null && "
      "qiime tools export --input-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+$x.qzv --output-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+$x >/dev/null"
      "; done ; done' > tmp/ANOSIM_{wildcards.sample}.sh &&"
      "chmod +x tmp/ANOSIM_{wildcards.sample}.sh &&"
      "bash tmp/ANOSIM_{wildcards.sample}.sh"




rule permdisp: # Calculates whether the two groups have similar dispersions (variances) to their centroid
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.beta_div.output.qza
   output:
                  temporary(expand("data/distance/PERMDISP/PERMDISP_{{sample}}+{group}+{dist}.qzv", dist=config["distances"], group=config["group"]))
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"])
   message: "Calculating PERMDISP for {wildcards.sample}"
   shell:
      "mkdir -p data/distance/PERMDISP &&"
      "echo 'for x in {params.dist}; do for y in {params.group}; do qiime diversity beta-group-significance --p-pairwise --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column $y --p-method permdisp --o-visualization data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+$x.qzv >/dev/null && "
      "qiime tools export --input-path data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+$x.qzv --output-path data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+$x >/dev/null"
      "; done ; done' > tmp/PERMDISP_{wildcards.sample}.sh &&"
      "chmod +x tmp/PERMDISP_{wildcards.sample}.sh &&"
      "bash tmp/PERMDISP_{wildcards.sample}.sh"



rule adonis: # Calculates whether the two groups have different centroids, susceptible to the groups having different dispersions. Therefore, interpret along with PERDISP (also known as betadisper)
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      rules.beta_div.output.qza
   output:
            temporary(expand("data/distance/ADONIS/ADONIS_{{sample}}+{group}+{dist}.qzv", dist=config["distances"], group=config["group"]))
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"])
   message: "Calculating ADONIS for {wildcards.sample}"
   shell:
      "mkdir -p data/distance/ADONIS &&"
      "echo 'for x in {params.dist}; do for y in {params.group}; do qiime diversity adonis --i-distance-matrix data/distance/beta_div/{wildcards.sample}+$x.qza --m-metadata-file data/map/{wildcards.sample}.txt --p-formula \"$y\" --o-visualization data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+$x.qzv >/dev/null && "
      "qiime tools export --input-path data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+$x.qzv --output-path data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+$x >/dev/null"
      "; done ; done' > tmp/ADONIS_{wildcards.sample}.sh &&"
      "chmod +x tmp/ADONIS_{wildcards.sample}.sh &&"
      "bash tmp/ADONIS_{wildcards.sample}.sh"


rule make_anosim_PDFs: # Creates PDFs out of the HTML file for easier handling
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.anosim.output
   output: 
      anosim=report(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+{dist}.pdf", dist=config["distances"], group=config["group"]))
   log: "data/logs/PDF_{sample}+ANOSIM.log"
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"])
   message: "Creating ANOSIM PDFs for {wildcards.sample}"
   shell:
      "mkdir -p data/logs &&"
      "echo 'for x in {params.dist}; do for y in {params.group}; do weasyprint data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+$x/index.html data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+$x.pdf 2>> data/logs/PDF_{wildcards.sample}+ANOSIM.log && rm -rf data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+$x"
      "; done ; done' > tmp/PDF_ANOSIM_{wildcards.sample}.sh &&"
      "chmod +x tmp/PDF_ANOSIM_{wildcards.sample}.sh &&"
      "bash tmp/PDF_ANOSIM_{wildcards.sample}.sh "


rule make_permdisp_PDFs: # Creates PDFs out of the HTML file for easier handling
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.permdisp.output
   output: 
      permdisp=report(expand("data/distance/PERMDISP/PERMDISP_{{sample}}+{group}+{dist}.pdf", dist=config["distances"], group=config["group"]))
   log: "data/logs/PDF_{sample}+PERMDISP.log"
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"])
   message: "Creating PERMDISP PDFs for {wildcards.sample}"
   shell:
      "mkdir -p data/logs &&"
      "echo 'for x in {params.dist}; do for y in {params.group}; do weasyprint data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+$x/index.html data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+$x.pdf 2>> data/logs/PDF_{wildcards.sample}+PERMDISP.log && rm -rf data/distance/PERMDISP/PERMDISP_{wildcards.sample}+$y+$x"
      "; done ; done' > tmp/PDF_PERMDISP_{wildcards.sample}.sh &&"
      "chmod +x tmp/PDF_PERMDISP_{wildcards.sample}.sh &&"
      "bash tmp/PDF_PERMDISP_{wildcards.sample}.sh "


rule make_adonis_PDFs: # Creates PDFs out of the HTML file for easier handling
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      rules.adonis.output
   output: 
      adonis=report(expand("data/distance/ADONIS/ADONIS_{{sample}}+{group}+{dist}.pdf", dist=config["distances"], group=config["group"]))
   log: "data/logs/PDF_{sample}+ADONIS.log"
   params: 
      dist=expand("{dist}",dist=config["distances"]),
      group=expand("{group}",group=config["group"])
   message: "Creating ADONIS PDFs for {wildcards.sample}"
   shell:
      "mkdir -p data/logs &&"
      "echo 'for x in {params.dist}; do for y in {params.group}; do weasyprint data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+$x/index.html data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+$x.pdf 2>> data/logs/PDF_{wildcards.sample}+ADONIS.log && rm -rf data/distance/ADONIS/ADONIS_{wildcards.sample}+$y+$x"
      "; done ; done' > tmp/PDF_ADONIS_{wildcards.sample}.sh &&"
      "chmod +x tmp/PDF_ADONIS_{wildcards.sample}.sh &&"
      "bash tmp/PDF_ADONIS_{wildcards.sample}.sh "
