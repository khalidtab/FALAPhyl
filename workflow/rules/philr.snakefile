rule correlation_clustering: #Opportunity for parallelization by writing all commands to a text file, then funnelling them to gnu parallel
   version: "1.0"
   conda:
      "../../workflow/envs/qiime2.yaml"
   input:
      "data/biom/{sample}.qza"
   output:
      folder=temporary(directory("data/tree/tree_{sample}")),
      tree_qza=temporary("data/tree/{sample}.qza"),
      tree=("data/tree/{sample}.tre")
   shell: 
      "mkdir -p data/tree data/plots &&"
      "qiime gneiss correlation-clustering --i-table {input} --o-clustering data/tree/{wildcards.sample}.qza &&"
      "qiime tools export --input-path data/tree/{wildcards.sample}.qza --output-path data/tree/tree_{wildcards.sample} &&"
      "mv data/tree/tree_{wildcards.sample}/tree.nwk data/tree/{wildcards.sample}.tre"

