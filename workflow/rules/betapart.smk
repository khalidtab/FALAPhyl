def get_mem_mb(wildcards, attempt):
    return attempt * 1000

rule betapart_matrix: # Create betapart matrices for all samples
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      repl=report("data/distance/beta_div/Repl–{sample}–{distance}.tsv",
      caption="../report/beta_matrix.rst",
      category="Beta diversity",
      subcategory="{distance}",
      labels={
         "Data type": "Distance Matrix",
         "Description": "Dissimilarity breakdown - Replacement component"}),
      norepl=report("data/distance/beta_div/NoRepl–{sample}–{distance}.tsv",
      caption="../report/beta_matrix.rst",
      category="Beta diversity",
      subcategory="{distance}",
      labels={
         "Data type": "Distance Matrix",
         "Description": "Dissimilarity breakdown - No replacement component",
         "Distance type": "{distance} breakdown"}),
      bray=temporary("data/distance/beta_div/{distance}_2–{sample}.tsv")
   message: "Creating the matrices for the {wildcards.distance} breakdown for {wildcards.sample}"
   shell:
      "Rscript --vanilla workflow/scripts/betapart_matrix.R -i {input} -r {output.repl} -n {output.norepl} -f {output.bray} -d {wildcards.distance}"

## 1. Pairwise analysis

checkpoint betapart_pairwise_write: # Create pairwise distance breakdown between each two entities in the categories to show whether the distance is mostly due to balanced variation in abundance or and balanced gradients
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      repl="data/distance/beta_div/Repl–{sample}–{distance}.tsv",
      norepl="data/distance/beta_div/NoRepl–{sample}–{distance}.tsv",
      bray="data/distance/beta_div/{distance}_2–{sample}.tsv"
   resources:
      mem_mb=get_mem_mb
   output:
      directory("data/plots/betapart_pairwise/{sample}–{distance}–{group}/")
   message: "Filtering {wildcards.distance} breakdown of pairwise distances for {wildcards.sample} {wildcards.group}"
   shell:
      "mkdir -p data/plots/betapart_pairwise/{wildcards.sample}–{wildcards.distance}–{wildcards.group}/ && "
      " Rscript --vanilla workflow/scripts/betapart_pairwise_write.R -i {input.bray} -r {input.repl} -n {input.norepl} -m data/{wildcards.sample}.txt -d {wildcards.distance} -c {wildcards.group} -o {output} "

rule betapart_pairwise_draw:
   conda:
      "../../workflow/envs/ggpubr.yaml"
   resources:
      mem_mb=get_mem_mb
   input:
      "data/plots/betapart_pairwise/{sample}–{distance}–{group}/{i}.tsv"
   output:
      report("data/plots/betapart_pairwise/{sample}–{distance}–{group}/{i}.svg",
      caption="../report/betapart_plot.rst",
      patterns=["{name}.svg"],
      category="Beta diversity pairwise",
      subcategory="{distance} breakdown",
      labels={
              "Description": "Betapart pairwise comparison",
              "Data type": "Plot",
              "Comparison": "{i}",
              "Method": "{distance} breakdown",
              "Grouping category": "{group}"})
   log: 
      "data/logs/betapart_pairwise_draw–{sample}–{distance}–{group}–{i}.txt"
   message: "Graphing {wildcards.distance} breakdown between {wildcards.i} for variable {wildcards.group} of {wildcards.sample}"
   shell:
      " Rscript --vanilla workflow/scripts/betapart_pairwise_draw.R -i {input} -m data/{wildcards.sample}.txt -d {wildcards.distance} -c {wildcards.group} -o {output} 2> {log} "

def pairwise_output(wildcards):
    checkpoint_output = checkpoints.betapart_pairwise_write.get(**wildcards).output[0]    
    return expand("data/plots/betapart_pairwise/{sample}–{distance}–{group}/{i}.svg",
           sample=wildcards.sample,
           group=wildcards.group,
           distance=wildcards.distance,
           i=glob_wildcards(os.path.join(checkpoint_output, "{i}.tsv")).i)


rule betapart_pairwise_done: # Create pairwise distance breakdown between each two entities in the categories to show whether the distance is mostly due to balanced variation in abundance or and balanced gradients
   input:
      pairwise_output
   output:
      touch("tmp/breakdown_pairwise_done_{sample}–{distance}–{group}.txt")
   message: "Done graphing breakdown of pairwise distances of {wildcards.distance} for {wildcards.sample} {wildcards.group}"

## 2. ANOSIM

rule repl_anosim_betapart:
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/distance/beta_div/Repl–{sample}–{distance}.tsv"
   log: 
      "data/logs/ANOSIM–Repl–{sample}–{distance}–{group}.txt"
   output:
      myreport=report("data/distance/ANOSIM/anosim–Repl–{sample}–{distance}–{group}.txt",
      caption="../report/beta_anosim.rst",
      category="Beta diversity",
      subcategory="{distance}",
      labels={
         "Data type": "Text file",
         "Description": "ANOSIM of the Replacement component",
         "Grouping category": "{group}",
         "Distance type": "{distance} breakdown"}),
      mypairwise=report("data/distance/ANOSIM/anosim–Repl–{sample}–{distance}–{group}_pairwise.txt",
      caption="../report/beta_anosim.rst",
      category="Beta diversity",
      subcategory="{distance}",
      labels={
         "Data type": "Text file",
         "Description": "ANOSIM of the Replacement component, pairwise comparisons with FDR correction",
         "Grouping category": "{group}",
         "Distance type": "{distance} breakdown"})
   message: "Calculating ANOSIM for {wildcards.distance} balanced (ie replacement) of {wildcards.sample} {wildcards.group}"
   shell:
      " Rscript --vanilla ./workflow/scripts/adonis_anosim_betadisper.R -i {input} -o {output.myreport} -m data/{wildcards.sample}.txt -g {wildcards.group} -t anosim > {log} "



use rule repl_anosim_betapart as norepl_anosim_betapart with:
   input:
      "data/distance/beta_div/NoRepl–{sample}–{distance}.tsv"
   log: 
      "data/logs/ANOSIM–NoRepl–{sample}–{distance}–{group}.txt"
   output:
      myreport=report("data/distance/ANOSIM/anosim–NoRepl–{sample}–{distance}–{group}.txt",
      caption="../report/beta_anosim.rst",
      category="Beta diversity",
      subcategory="{distance}",
      labels={
         "Data type": "Text file",
         "Description": "ANOSIM of the No replacement component",
         "Grouping category": "{group}",
         "Distance type": "{distance} breakdown"}),
      mypairwise=report("data/distance/ANOSIM/anosim–NoRepl–{sample}–{distance}–{group}_pairwise.txt",
      caption="../report/beta_anosim.rst",
      category="Beta diversity",
      subcategory="{distance}",
      labels={
         "Data type": "Text file",
         "Description": "ANOSIM of the No replacement component, pairwise comparisons with FDR correction",
         "Grouping category": "{group}",
         "Distance type": "{distance} breakdown"})         
   message: "Calculating ANOSIM for {wildcards.distance} gradient (ie no replacement) of {wildcards.sample} {wildcards.group}"

## 3. Graphing breakdown based on permutation

### Split biom file for the upcoming betapart steps

checkpoint group_and_category: # Will export multiple biom files tagged by their group name (retrieved from parameter) and the category within that group ("group+category.txt"). The tagging is done within the R script
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory("data/tsv/{sample}–{group}/tsv/"))
   message: "Creating biom files separated by {wildcards.group} for {wildcards.sample}"
   shell:
      "mkdir -p data/tsv/{wildcards.sample}–{wildcards.group}/tsv/ && "
      "Rscript --vanilla workflow/scripts/separate_bioms.R -i {input} -o {output} -m data/{wildcards.sample}.txt -g {wildcards.group} "

### 3.A
rule ingroup_var: # Create betapart mean, standard deviation
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}–{group}/tsv/{id}.tsv"
   output:
      ingroupdiff=report("data/betapart/{sample}–{distance}–{group}/ingroupdiff/{id}–ingroupdiff.txt",
      category="Beta diversity",
      subcategory="{distance}",
            labels={
         "Data type": "Text file",
         "Description": "In-group variance calc. – {id}",
         "Grouping category": "{group}",
         "Distance type": "{distance} breakdown"}),
      permu_file="data/betapart/{sample}–{distance}–{group}/permutations/{id}–permutations.txt"
   message: "Calculating mean and standard deviation, and permutation based distribution for {wildcards.distance} within {wildcards.id} {wildcards.group}."
   params:
      per_perm=expand("{per_perm}",per_perm=config["betapart_samples"]),
      permutations=expand("{permutations}",permutations=config["betapart_permutations"])
   log: 
      "data/logs/ingroupdiff–{sample}–{distance}–{group}–{id}.log"
   shell:
      " Rscript --vanilla workflow/scripts/betapart_per_id.R -i {input} -z {output.ingroupdiff} -d {wildcards.distance} -r {params.per_perm} -p {params.permutations} -o {output.permu_file} > {log} "

def ids_ingroup_var(wildcards):
    splitbiom_checkpoint_output = checkpoints.group_and_category.get(**wildcards).output[0]    
    file_names = expand("data/betapart/{sample}–{distance}–{group}/ingroupdiff/{id}–ingroupdiff.txt", sample=wildcards.sample, group=wildcards.group, distance=wildcards.distance, id = glob_wildcards(os.path.join(splitbiom_checkpoint_output, "{id}.tsv")).id)
    return file_names
                         
rule ngroup_var_done: # Mean standard deviation of breakdown is done
   input:
      ids_ingroup_var
   output:
      touch("tmp/{sample}–{distance}–{group}–ingroupdiff.txt")
   message: "Done with calculating mean and standard deviation for {wildcards.distance}."


### 3.B
rule betapart_plot_permutation: # Plots permutations from betapart_permutations step
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      ids_ingroup_var
   output:
      report("data/plots/Breakdown–{sample}–{distance}–{group}.svg",
      category="Beta diversity",
      subcategory="{distance}",
            labels={
         "Data type": "Plot - SVG",
         "Description": "Permutation-based dissimilarity breakdown",
         "Grouping category": "{group}",
         "Distance type": "{distance} breakdown"})
   params:
      color=expand("{color}",color=config["color"]),
      width=config["width"],
      height=config["height"]
   message: "Plotting permutations for {wildcards.distance} breakdown of {wildcards.sample} {wildcards.group}"
   log: "data/logs/betapart_permutation_plot_{sample}–{distance}–{group}.log"
   shell:
      " Rscript --vanilla workflow/scripts/betapart_plot_permutation.R -i data/betapart/{wildcards.sample}–{wildcards.distance}–{wildcards.group}/permutations/ -o {output} -m data/{wildcards.sample}.txt -b {wildcards.sample} -c {params.color} -d {wildcards.distance} -x {params.width} -y {params.height} > {log} "


### Finalizing betapart results

rule betapart:
   input:
      "data/distance/ANOSIM/anosim–Repl–{sample}–{distance}–{group}.txt",
      "data/distance/ANOSIM/anosim–NoRepl–{sample}–{distance}–{group}.txt",
      "tmp/breakdown_pairwise_done_{sample}–{distance}–{group}.txt",
      "data/plots/Breakdown–{sample}–{distance}–{group}.svg",
      "tmp/{sample}–{distance}–{group}–ingroupdiff.txt"
   output:
      touch("tmp/betapart–{distance}–{sample}–{group}.txt")


rule breakdown:
   input:
      expand("tmp/betapart–{distance}–{sample}–{group}.txt",  sample=config["mysample"],group=config["group"],distance=["bray","jaccard"])