rule betapart_matrix: # Create betapart matrices for all samples
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      "data/biom/{sample}.biom"
   output:
      repl="data/distance/beta_div/{sample}+BrayRepl.tsv",
      norepl="data/distance/beta_div/{sample}+BrayNoRepl.tsv"
   shell:
      "mkdir -p data/distance/beta_div/{wildcards.sample} && Rscript --vanilla workflow/scripts/betapart_matrix.R -i {input} -r {output.repl} -n {output.norepl}"



rule anosim_betapart: # Calculates whether the intra-group variances is sig different from intergroup variances
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      repl  =rules.betapart_matrix.output.repl,
      norepl=rules.betapart_matrix.output.norepl
   output:
      BrayReplDist=temporary("data/distance/beta_div/{{sample}}+BrayRepl.qza"),
      BrayNoReplDist=temporary("data/distance/beta_div/{{sample}}+BrayNoRepl.qza"),
      BrayRepl=temporary(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+BrayRepl.qzv", group=config["group"])),
      BrayNoRepl=temporary(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+BrayNoRepl.qzv", group=config["group"]))
   params: 
      group=expand("{group}",group=config["group"])
   message: "Calculating ANOSIM for betapart of {wildcards.sample}"
   shell:
      "mkdir -p data/distance/ANOSIM && mkdir -p data/distance/ANOSIM && "
      "qiime tools import --input-path {input.repl} --output-path {output.BrayReplDist} --type DistanceMatrix &&"
      "qiime tools import --input-path {input.norepl} --output-path {output.BrayNoReplDist} --type DistanceMatrix &&"
      "echo 'for y in {params.group}; do "
      "qiime diversity beta-group-significance --p-pairwise --i-distance-matrix {output.BrayReplDist} --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column $y --p-method anosim --o-visualization data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayRepl.qzv >/dev/null && "
      "qiime diversity beta-group-significance --p-pairwise --i-distance-matrix {output.BrayNoReplDist} --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column $y --p-method anosim --o-visualization data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayNoRepl.qzv >/dev/null && "
      "qiime tools export --input-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayRepl.qzv --output-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayRepl >/dev/null &&"
      "qiime tools export --input-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayNoRepl.qzv --output-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayNoRepl >/dev/null"
      "; done' > tmp/ANOSIM_Betapart_{wildcards.sample}.sh &&"
      "chmod +x tmp/ANOSIM_Betapart_{wildcards.sample}.sh &&"
      "bash tmp/ANOSIM_Betapart_{wildcards.sample}.sh"


rule make_anosim_betapart_PDFs: # Create betapart matrices for all samples
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      BrayRepl=rules.anosim_betapart.output.BrayRepl,
      BrayNoRepl=rules.anosim_betapart.output.BrayNoRepl
   output: 
      anosim_norepl=report(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+BrayNoRepl.pdf", group=config["group"])),
      anosim_repl=report(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+BrayRepl.pdf", group=config["group"]))
   log: "data/logs/PDF_{sample}+ANOSIM+Betapart.log"
   params: 
      group=expand("{group}",group=config["group"])
   message: "Creating ANOSIM PDFs for Betapart of {wildcards.sample}"
   shell:
      "mkdir -p data/logs &&"
      "echo 'for y in {params.group}; do weasyprint data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayNoRepl/index.html data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayNoRepl.pdf 2>> data/logs/PDF_{wildcards.sample}+ANOSIM+Betapart.log && rm -rf data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayNoRepl && "
      "weasyprint data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayRepl/index.html data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayRepl.pdf 2>> data/logs/PDF_{wildcards.sample}+ANOSIM+Betapart.log && rm -rf data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+BrayRepl"
      "; done' > tmp/PDF_ANOSIM_betapart_{wildcards.sample}.sh &&"
      "chmod +x tmp/PDF_ANOSIM_betapart_{wildcards.sample}.sh &&"
      "bash tmp/PDF_ANOSIM_betapart_{wildcards.sample}.sh "

# Split biom file for the upcoming betapart steps

checkpoint group_and_category: # Will export multiple biom files tagged by their group name (retrieved from parameter) and the category within that group ("group+category.txt"). The tagging is done within the R script
   version: "1.0"
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      "data/biom/{sample}.biom"
   params:
      group=expand("{group}",group=config["group"])
   output:
      temporary(directory("data/betapart/{sample}/tsv/"))
   message: "Creating biom files separated by {params.group} for {wildcards.sample}"
   shell:
      "mkdir -p data/betapart/{wildcards.sample}/tsv data/logs && "
      "echo 'for x in {params.group}; do "
      "Rscript --vanilla workflow/scripts/separate_bioms.R -i data/biom/{wildcards.sample}.biom -o data/betapart/{wildcards.sample}/tsv -m data/map/{wildcards.sample}.txt -g $x ; done ' > tmp/getbiomsPerCategory_{wildcards.sample}.sh &&"
      "chmod +x tmp/getbiomsPerCategory_{wildcards.sample}.sh &&"
      "bash tmp/getbiomsPerCategory_{wildcards.sample}.sh"



rule ingroup_var: # Create betapart mean, standard deviation
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      "data/betapart/{sample}/tsv/{id}.tsv"
   output:
      "data/betapart/{sample}/ingroupdiff/{id}_ingroupdiff.txt"
   message: "Calculating mean and standard deviation within {wildcards.id}."
   shell:
      "mkdir -p data/betapart/{wildcards.sample}/ingroupdiff && touch {output} && Rscript --vanilla workflow/scripts/betapart_ingroup_var.R -i {input} -o {output}"
      

def ids_ingroup_var(wildcards):
    splitbiom_checkpoint_output = checkpoints.group_and_category.get(**wildcards).output[0]    
    file_names = expand("data/betapart/{sample}/ingroupdiff/{id}_ingroupdiff.txt", sample=wildcards.sample, id = glob_wildcards(os.path.join(splitbiom_checkpoint_output, "{id}.tsv")).id)
    return file_names





rule betapart_permutations: # Create distribution per permutation
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      "data/betapart/{sample}/tsv/{id}.tsv"
   output:
      temporary(touch("data/betapart/{sample}/perm/{id}_betapart_permutation.txt"))
   params:
      per_perm=expand("{per_perm}",per_perm=config["betapart_samples"]),
      permutations=expand("{permutations}",permutations=config["betapart_permutations"])
   log: "data/logs/betapart_permutations+{sample}+{id}.log"
   message: "Starting permutations for {wildcards.id}…"
   shell:
      "mkdir -p data/betapart/{wildcards.sample}/perm/permutations data/logs/ && Rscript --vanilla workflow/scripts/betapart_permutation.R -i {input} -r {params.per_perm} -p {params.permutations} -o data/betapart/{wildcards.sample}/perm >> data/logs/betapart_permutations+{wildcards.sample}+{wildcards.id}.log"
      


def ids_betapart_permutations(wildcards):
    splitbiom_checkpoint_output = checkpoints.group_and_category.get(**wildcards).output[0]    
    file_names = expand("data/betapart/{sample}/perm/{id}_betapart_permutation.txt", sample=wildcards.sample, id = glob_wildcards(os.path.join(splitbiom_checkpoint_output, "{id}.tsv")).id)
    return file_names

rule betapart_plot_permutation: # Plots permutations from betapart_permutations step
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      perm=ids_betapart_permutations
   output:
      temporary(touch("data/betapart/{sample}/plots_done.txt"))
   params:
      color=expand("{color}",color=config["color"])
   log: "data/logs/PDF_{sample}+ANOSIM+Betapart.log"
   shell:
      "mkdir -p data/plots/ && Rscript --vanilla workflow/scripts/betapart_plot_permutation.R -i data/betapart/{wildcards.sample}/perm/permutations/ -o data/plots/ -m data/map/{wildcards.sample}.txt -b {wildcards.sample} -c {params.color}"





rule betapart: # Last step from betapart. Cleans up temporary files.
    input:
        mean_std=ids_ingroup_var,
        perm=rules.betapart_plot_permutation.output,
        PDFs=rules.make_anosim_betapart_PDFs.output.anosim_repl
    output:
        temporary(touch("data/.done_ingroupdiff_{sample}.txt"))
    shell:
        "echo Cleaning up after Betapart. && rm -rf data/betapart/{wildcards.sample}/perm/permutations data/betapart/{wildcards.sample}/tsv data/distance/beta_div/{wildcards.sample}"
