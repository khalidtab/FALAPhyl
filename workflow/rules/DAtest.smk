checkpoint biom_pairwise: # Create pairwise biom files for all the variables so that the differential abundance tests would work
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory("data/diff/pairwise–{sample}–{group}/"))
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}"
   shell:
      " mkdir -p {output} && "
      " Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o {output} "

rule testDA:
   version: "1.0"
   container:
        "docker://khalidtab/datest:latest"
   input:
      tsv="data/diff/pairwise–{sample}–{group}/{annie}.tsv",
      map="data/diff/pairwise–{sample}–{group}/{annie}.txt"
   output:
      testfile="data/diff_{sample}–{group}/EmPower_FDR_{annie}.tsv",
      graph="data/diff_{sample}–{group}/EmPower_FDR_{annie}.svg"
   message: "testDA started {wildcards.sample} – {wildcards.group} for {wildcards.annie}"
   shell:
      "Rscript --vanilla workflow/scripts/testDA.R -i {input.tsv} -m {input.map} -c {wildcards.group} -o {output.testfile} -g {output.graph}"


def ids_testDA(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise.get(**wildcards).output[0]    
    file_names = expand("data/diff_{sample}–{group}/EmPower_FDR_{annie}.svg", sample=wildcards.sample, group=wildcards.group, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names
                         
rule testDA_done: 
   input:
      ids_testDA
   output:
      touch("tmp/testDA_{sample}–{group}–done.txt")


rule diff:
   input:
      expand("tmp/testDA_{sample}–{group}–done.txt",  sample=config["mysample"],group=config["group"])