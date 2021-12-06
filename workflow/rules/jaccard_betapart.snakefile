rule jaccard_betapart_matrix: # Create betapart matrices for all samples
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      repl="data/distance/beta_div/{sample}+jaccardRepl.tsv",
      norepl="data/distance/beta_div/{sample}+jaccardNoRepl.tsv",
      jaccard=temporary("data/distance/beta_div/{sample}+jaccard_2.tsv")
   message: "Creating the matrices for the Jaccard breakdown for {wildcards.sample}"
   shell:
      "mkdir -p data/distance/beta_div/{wildcards.sample} && Rscript --vanilla workflow/scripts/betapart_matrix.R -i {input} -r {output.repl} -n {output.norepl} -f {output.jaccard} -d jaccard"

rule jaccard_betapart_pairwise: # Create pairwise distance breakdown between each two entities in the categories to show whether the distance is mostly due to turn-over or nestedness
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      repl="data/distance/beta_div/{sample}+jaccardRepl.tsv",
      norepl="data/distance/beta_div/{sample}+jaccardNoRepl.tsv",
      jaccard="data/distance/beta_div/{sample}+jaccard_2.tsv"
   params:
      group=expand("{group}", group=config["group"])
   output:
      temporary(touch("data/distance/beta_div/{sample}_jaccard_breakdown_pairwise.done"))
   message: "Graphing breakdown of pairwise distances of Jaccard for {wildcards.sample}"
   shell:
      "mkdir -p data/distance/beta_div/{wildcards.sample} &&"
      "echo 'for y in {params.group}; do "
      " Rscript --vanilla workflow/scripts/betapart_pairwise.R -i {input.jaccard} -r {input.repl} -n {input.norepl} -m data/map/{wildcards.sample}.txt -d jaccard -c $y -o data/plots/betapart_jaccard_pairwise_{wildcards.sample}_"
      " ; done' > tmp/Betapart_jaccard_pairwise_{wildcards.sample}.sh &&"
      " chmod +x tmp/Betapart_jaccard_pairwise_{wildcards.sample}.sh &&"
      " bash tmp/Betapart_jaccard_pairwise_{wildcards.sample}.sh"

rule jaccard_anosim_betapart: # Calculates whether the intra-group variances is sig different from intergroup variances
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      repl  =rules.jaccard_betapart_matrix.output.repl,
      norepl=rules.jaccard_betapart_matrix.output.norepl
   output:
      jaccardReplDist=temporary("data/distance/beta_div/{{sample}}+jaccardRepl.qza"),
      jaccardNoReplDist=temporary("data/distance/beta_div/{{sample}}+jaccardNoRepl.qza"),
      jaccardRepl=temporary(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+jaccardRepl.qzv", group=config["group"])),
      jaccardNoRepl=temporary(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+jaccardNoRepl.qzv", group=config["group"]))
   params: 
      group=expand("{group}",group=config["group"])
   message: "Calculating ANOSIM for Jaccard betapart of {wildcards.sample}"
   shell:
      "mkdir -p data/distance/ANOSIM && mkdir -p data/distance/ANOSIM && "
      "qiime tools import --input-path {input.repl} --output-path {output.jaccardReplDist} --type DistanceMatrix &&"
      "qiime tools import --input-path {input.norepl} --output-path {output.jaccardNoReplDist} --type DistanceMatrix &&"
      "echo 'for y in {params.group}; do "
      "qiime diversity beta-group-significance --p-pairwise --i-distance-matrix {output.jaccardReplDist} --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column $y --p-method anosim --o-visualization data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardRepl.qzv >/dev/null && "
      "qiime diversity beta-group-significance --p-pairwise --i-distance-matrix {output.jaccardNoReplDist} --m-metadata-file data/map/{wildcards.sample}.txt --m-metadata-column $y --p-method anosim --o-visualization data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardNoRepl.qzv >/dev/null && "
      "qiime tools export --input-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardRepl.qzv --output-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardRepl >/dev/null &&"
      "qiime tools export --input-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardNoRepl.qzv --output-path data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardNoRepl >/dev/null"
      "; done' > tmp/ANOSIM_Betapart_jaccard_{wildcards.sample}.sh &&"
      "chmod +x tmp/ANOSIM_Betapart_jaccard_{wildcards.sample}.sh &&"
      "bash tmp/ANOSIM_Betapart_jaccard_{wildcards.sample}.sh"


rule jaccard_make_anosim_betapart_PDFs: # Create betapart matrices for all samples
   version: "1.0"
   conda: "../../workflow/envs/html2pdf.yaml"
   input:
      jaccardRepl=rules.jaccard_anosim_betapart.output.jaccardRepl,
      jaccardNoRepl=rules.jaccard_anosim_betapart.output.jaccardNoRepl
   output: 
      anosim_norepl=report(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+jaccardNoRepl.pdf", group=config["group"])),
      anosim_repl=report(expand("data/distance/ANOSIM/ANOSIM_{{sample}}+{group}+jaccardRepl.pdf", group=config["group"]))
   log: "data/logs/PDF_{sample}+ANOSIM+jaccard_Betapart.log"
   params: 
      group=expand("{group}",group=config["group"])
   message: "Creating ANOSIM PDFs for Jaccard betapart of {wildcards.sample}"
   shell:
      "mkdir -p data/logs &&"
      "echo 'for y in {params.group}; do weasyprint data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardNoRepl/index.html data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardNoRepl.pdf 2>> data/logs/PDF_{wildcards.sample}+ANOSIM+Betapart.log && rm -rf data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardNoRepl && "
      "weasyprint data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardRepl/index.html data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardRepl.pdf 2>> data/logs/PDF_{wildcards.sample}+ANOSIM+Betapart.log && rm -rf data/distance/ANOSIM/ANOSIM_{wildcards.sample}+$y+jaccardRepl"
      "; done' > tmp/PDF_ANOSIM_jaccard_betapart_{wildcards.sample}.sh &&"
      "chmod +x tmp/PDF_ANOSIM_jaccard_betapart_{wildcards.sample}.sh &&"
      "bash tmp/PDF_ANOSIM_jaccard_betapart_{wildcards.sample}.sh "

# Split biom file for the upcoming betapart steps

checkpoint jaccard_group_and_category: # Will export multiple biom files tagged by their group name (retrieved from parameter) and the category within that group ("group+category.txt"). The tagging is done within the R script
   version: "1.0"
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      "data/tsv/{sample}.tsv"
   params:
      group=expand("{group}",group=config["group"])
   output:
      temporary(directory("data/betapart_jaccard/{sample}/tsv/"))
   message: "Creating biom files separated by {params.group} for {wildcards.sample} for Jaccard"
   shell:
      "mkdir -p data/betapart_jaccard/{wildcards.sample}/tsv data/logs && "
      "echo 'for x in {params.group}; do "
      "Rscript --vanilla workflow/scripts/separate_bioms.R -i data/tsv/{wildcards.sample}.tsv -o data/betapart_jaccard/{wildcards.sample}/tsv -m data/map/{wildcards.sample}.txt -g $x ; done ' > tmp/getbiomsPerCategoryjaccard_{wildcards.sample}.sh &&"
      "chmod +x tmp/getbiomsPerCategoryjaccard_{wildcards.sample}.sh &&"
      "bash tmp/getbiomsPerCategoryjaccard_{wildcards.sample}.sh"



rule jaccard_ingroup_var: # Create betapart mean, standard deviation
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      "data/betapart_jaccard/{sample}/tsv/{id}.tsv"
   output:
      "data/betapart_jaccard/{sample}/ingroupdiff/{id}_ingroupdiff.txt"
   message: "Calculating mean and standard deviation for Jaccard within {wildcards.id}."
   shell:
      "mkdir -p data/betapart_jaccard/{wildcards.sample}/ingroupdiff && touch {output} && Rscript --vanilla workflow/scripts/betapart_ingroup_var.R -i {input} -o {output} -d jaccard"
      

def jaccard_ids_ingroup_var(wildcards):
    splitbiom_checkpoint_output = checkpoints.jaccard_group_and_category.get(**wildcards).output[0]    
    file_names = expand("data/betapart_jaccard/{sample}/ingroupdiff/{id}_ingroupdiff.txt", sample=wildcards.sample, id = glob_wildcards(os.path.join(splitbiom_checkpoint_output, "{id}.tsv")).id)
    return file_names





rule betapart_jaccard_permutations: # Create distribution per permutation
   conda:
      "../../workflow/envs/betapart.yaml"
   input:
      "data/betapart_jaccard/{sample}/tsv/{id}.tsv"
   output:
      temporary(touch("data/betapart_jaccard/{sample}/perm/{id}_betapart_permutation.txt"))
   params:
      per_perm=expand("{per_perm}",per_perm=config["betapart_samples"]),
      permutations=expand("{permutations}",permutations=config["betapart_permutations"])
   log: "data/logs/betapart_permutations+{sample}+{id}.log"
   message: "Starting Jaccard permutations for {wildcards.id}…"
   shell:
      "mkdir -p data/betapart_jaccard/{wildcards.sample}/perm/permutations data/logs/ && Rscript --vanilla workflow/scripts/betapart_permutation.R -i {input} -r {params.per_perm} -p {params.permutations} -d jaccard -o data/betapart_jaccard/{wildcards.sample}/perm/ >> data/logs/betapart_permutations+{wildcards.sample}+{wildcards.id}.log"
      


def ids_jaccard_betapart_permutations(wildcards):
    splitbiom_checkpoint_output = checkpoints.jaccard_group_and_category.get(**wildcards).output[0]    
    file_names = expand("data/betapart_jaccard/{sample}/perm/{id}_betapart_permutation.txt", sample=wildcards.sample, id = glob_wildcards(os.path.join(splitbiom_checkpoint_output, "{id}.tsv")).id)
    return file_names

rule jaccard_betapart_plot_permutation: # Plots permutations from betapart_permutations step
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      perm=ids_jaccard_betapart_permutations
   output:
      temporary(touch("data/betapart_jaccard/{sample}/plots_done.txt"))
   params:
      color=expand("{color}",color=config["color"])
   message: "Creating the permutations for the Jaccard breakdown of {wildcards.sample}"
   log: "data/logs/PDF_{sample}+ANOSIM+Betapart.log"
   shell:
      "mkdir -p data/plots/ && Rscript --vanilla workflow/scripts/betapart_plot_permutation.R -i data/betapart_jaccard/{wildcards.sample}/perm/permutations/ -o data/plots/jaccard_ -m data/map/{wildcards.sample}.txt -b {wildcards.sample} -c {params.color} -d jaccard"





rule jaccard_betapart: # Last step from betapart. Cleans up temporary files.
    input:
        mean_std=jaccard_ids_ingroup_var,
        perm=rules.jaccard_betapart_plot_permutation.output,
        PDFs=rules.jaccard_make_anosim_betapart_PDFs.output.anosim_repl,
        pairwise=rules.jaccard_betapart_pairwise.output
    output:
        temporary(touch("data/.done_jaccard_betapart_{sample}.txt"))
    shell:
        "echo Cleaning up after Jaccard Betapart. && rm -rf data/betapart_jaccard/{wildcards.sample}/perm/permutations data/betapart_jaccard/{wildcards.sample}/tsv data/distance/beta_div/{wildcards.sample}"


