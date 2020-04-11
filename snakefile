configfile: "input.yaml"
#container: "docker://continuumio/miniconda3:4.8.2"
include: "workflow/rules/make_qza.snakefile"
include: "workflow/rules/beta_div.snakefile"

#input: lambda wildcards: glob('sample-{samp}/*.fastq'.format(samp=wildcards.samp))

rule all:
   input: "{sample}.final"
   shell: "touch {input}"


rule create_stop_trigger:
   version: "1.0"
   conda: "workflow/envs/qiime2.yaml"
   input:
      rules.make_anosim_PDFs.output,
      rules.make_permdisp_PDFs.output,
      rules.make_adonis_PDFs.output,
      rules.pcoa_svg.output
   output: 
      "{sample}.final"
   shell:
      "touch {output}"

