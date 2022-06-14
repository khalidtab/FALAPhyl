rule biom_to_tsv: # Since qiime requires a specific biom format (either V100/json or V210/hdf5), we can take the generic biom file, and force it to be a json, then import it to qiime as an artifact
   version: "1.0"
   conda:
        "../../workflow/envs/biom.yaml"
   input:
        "data/{sample}.biom"
   output:
        tsv=temporary("data/tsv/{sample}.tsv"),
        temp_biom=temporary("data/biom/{sample}_temp.biom")
   message: "Generating temporary tsv files of the biom file for {wildcards.sample}"
   shell:
        "mkdir -p data/tsv &&"
        "biom convert -i {input} -o {output.tsv} --to-tsv &&"
        'biom convert -i {output.tsv} -o {output.temp_biom} --to-json --table-type="OTU table"'