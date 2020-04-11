rule biom_to_qza: # Since qiime requires a specific biom format (either V100/json or V210/hdf5), we can take the generic biom file, and force it to be a json, then import it to qiime as an artifact
   version: "1.0"
   conda: "../../workflow/envs/qiime2.yaml"
   input:
      "data/biom/{sample}.biom"
   output:
        qza=temporary("data/biom/{sample}.qza"),
        tsv=temporary("data/tsv/{sample}.tsv"),
        temp_biom=temporary("data/biom/{sample}_temp.biom")
   shell:
      "mkdir -p data/tsv &&"
      "biom convert -i {input} -o {output.tsv} --to-tsv --header-key taxonomy &&"
      'biom convert -i {output.tsv} -o {output.temp_biom} --to-json --table-type="OTU table" --process-obs-metadata taxonomy &&'
      "qiime tools import --input-path {output.temp_biom} --output-path {output.qza} --type FeatureTable[Frequency] --input-format BIOMV100Format"

