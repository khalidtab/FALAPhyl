rule biom_to_tsv: # Since qiime requires a specific biom format (either V100/json or V210/hdf5), we can take the generic biom file, and force it to be a json, then import it to qiime as an artifact
   version: "1.0"
   conda:
        "../../workflow/envs/biom.yaml"
   input:
        "data/{sample}.biom"
   output:
        "data/tsv/{sample}_full.tsv"
   message: "Generating temporary tsv files of the biom file for {wildcards.sample}"
   shell:
      "mkdir -p data/tsv &&"
      "biom convert -i {input} -o {output} --to-tsv "

rule filter_biom: # Remove samples that are not in the mapping file
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}_full.tsv"
   output:
      "data/tsv/{sample}.tsv"
   message: "Removing any samples from the biom file that is not part of the mapping file for {wildcards.sample}"
   shell: 
      "Rscript --vanilla ./workflow/scripts/filter_biom.R -i {input} -o {output} -m data/{wildcards.sample}.txt "

rule create_filtered_biom: 
   version: "2.0"
   conda:
      "../../workflow/envs/biom.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      "data/biom/{sample}_temp.biom"
   log:
      "data/logs/filterBiom–{sample}.txt"
   message: "Creating filtered biom file for {wildcards.sample}"
   shell: 
      'biom convert -i {input} -o {output} --to-json --table-type="OTU table"  > {log} '


checkpoint biom_pairwise1: # for DA effect size calculation
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory(expand("data/tmp/diff/{{sample}}–{{group}}–minAbd{{minabund}}minR{{minread}}minS{{minsample}}/{theTest}–pairwise/",theTest=config["DA_tests"])))
   params: theTest=config['DA_tests']
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}"
   shell:
      '''
      mkdir -p {output} data/tmp/difftmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}
      Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o data/tmp/difftmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}
      # Set IFS to newline to handle file paths with special characters or spaces
      IFS=$' '
      # List of unquoted items (simulating unquoted terminal output)
      items=({params.theTest})
      
      # Loop through each item in the list
      for item in ${{items[@]}}; do
          cp -r data/tmp/difftmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/* data/tmp/diff/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/$item–pairwise/
      done
      rm -r data/tmp/difftmp/
      # Reset IFS to default (optional)
      unset IFS
      '''

checkpoint biom_pairwise2: # DAtest
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory("data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/"))
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}"
   shell:
      " mkdir -p {output} | "
      " Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o {output} "

