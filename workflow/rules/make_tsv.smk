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

rule include_biom_and_meta:
   conda:
      "../../workflow/envs/biom.yaml"
   input:
      expand("data/tsv/{sample}.tsv",sample=config["mysample"])
   output:
      biom=report(expand("data/biom/{sample}_temp.biom",sample=config["mysample"]),
             category="Biom file",
             subcategory="Input files",
             labels={"Data type": "Count table in BIOM-format"}),
      metadata=report(expand("data/tmp/{sample}.txt",sample=config["mysample"]),
             category="Mapping file",
             subcategory="Input files",
             labels={"Data type": "Mapping file"})
   params: mySample=config["mysample"]
   shell:
      'biom convert -i {input} -o {output.biom} --to-json --table-type="OTU table" && cp data/{params.mySample}.txt data/tmp/{params.mySample}.txt'

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

