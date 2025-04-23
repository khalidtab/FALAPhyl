rule biom_to_tsv: # Since qiime requires a specific biom format (either V100/json or V210/hdf5), we can take the generic biom file, and force it to be a json, then import it to qiime as an artifact
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


checkpoint makeBiomForEffectSize: # EffectSize
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory(expand("data/tmp/diff/EffectSize–{{sample}}–{{group}}–minAbd{{minabund}}minR{{minread}}minS{{minsample}}–{theTest}/",theTest=config["DA_tests"])))
   params: theTest=config['DA_tests']
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}, and creating temporary files for each test to avoid any file access clashes"
   shell:
      '''
      mkdir -p {output} data/tmp/difftmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}
      Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o data/tmp/difftmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}
      IFS=$"\n"  # Newline-separated list, in case test names contain spaces
      
      items=({params.theTest})
      
      # Loop through each DA test method
      for item in "${{items[@]}}"; do
      # Create target directory for the test
      target_dir="data/tmp/diff/EffectSize–{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}–$item"
      mkdir -p "$target_dir"
      
      # Loop through each file in the temp diff directory
      for filepath in data/tmp/difftmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/*; do
      filename_with_ext=$(basename "$filepath")
      extension="${{filename_with_ext##*.}}"      # File extension
      filename="${{filename_with_ext%.*}}"        # Filename without extension
      
      # Construct new name and copy
      cp "$filepath" "$target_dir/${{filename}}.${{extension}}"
      done
      done
      
      # Cleanup
      rm -r data/tmp/difftmp/
      unset IFS
      '''

checkpoint makeBiomForEffectSizePaired: # EffectSize
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory(expand("data/tmp/diffPaired/EffectSize–{{sample}}–{{group}}–minAbd{{minabund}}minR{{minread}}minS{{minsample}}–{theTest}/",theTest=config["paired_DA_tests"])))
   params: theTest=config['paired_DA_tests']
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}, and creating temporary files for each test to avoid any file access clashes"
   shell:
      '''
      mkdir -p {output} data/tmp/diffPairedtmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}
      Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o data/tmp/diffPairedtmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}
      IFS=$"\n"  # Newline-separated list, in case test names contain spaces
      
      items=({params.theTest})
      
      # Loop through each DA test method
      for item in "${{items[@]}}"; do
      # Create target directory for the test
      target_dir="data/tmp/diffPaired/EffectSize–{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}–$item"
      mkdir -p "$target_dir"
      
      # Loop through each file in the temp diff directory
      for filepath in data/tmp/diffPairedtmp/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/*; do
      filename_with_ext=$(basename "$filepath")
      extension="${{filename_with_ext##*.}}"      # File extension
      filename="${{filename_with_ext%.*}}"        # Filename without extension
      
      # Construct new name and copy
      cp "$filepath" "$target_dir/${{filename}}.${{extension}}"
      done
      done
      
      # Cleanup
      rm -r data/tmp/diffPairedtmp/
      unset IFS
      '''

checkpoint makeBiomForDAtest: # DAtest
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory("data/tmp/diff/DAtest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/"))
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}"
   shell:
      " mkdir -p {output} | "
      " Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o {output} "

