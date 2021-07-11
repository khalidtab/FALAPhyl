rule biom_to_qza: # Since qiime requires a specific biom format (either V100/json or V210/hdf5), we can take the generic biom file, and force it to be a json, then import it to qiime as an artifact
   version: "1.0"
   conda:
        "../../workflow/envs/qiime2.yaml"
   input:
        "data/biom/{sample}.biom"
   output:
        qza=("data/biom/{sample}.qza"),
        tsv=("data/tsv/{sample}.tsv"),
        temp_biom=("data/biom/{sample}_temp.biom")
   message: "Generating temporary files for {wildcards.sample}"
   shell:
        "mkdir -p data/tsv &&"
        "biom convert -i {input} -o {output.tsv} --to-tsv &&"
        'biom convert -i {output.tsv} -o {output.temp_biom} --to-json --table-type="OTU table" &&'
        "qiime tools import --input-path {output.temp_biom} --output-path {output.qza} --type FeatureTable[Frequency] --input-format BIOMV100Format >/dev/null"