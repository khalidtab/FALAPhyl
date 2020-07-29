configfile: "data/input.yaml"
include: "workflow/rules/fastspar.snakefile"
include: "workflow/rules/make_qza.snakefile"
include: "workflow/rules/beta_div.snakefile"
include: "workflow/rules/alpha_div.snakefile"
include: "workflow/rules/beta_sig.snakefile"
include: "workflow/rules/philr.snakefile"

rule all:
   input: "{sample}.final"
   shell: "touch {input}"

rule create_stop_trigger:
   version: "1.0"
   conda: "workflow/envs/qiime2.yaml"
   input:
      rules.run_fastspar.output,
      rules.make_anosim_PDFs.output,
      rules.make_permdisp_PDFs.output,
      rules.make_adonis_PDFs.output,
      rules.pcoa_svg.output,
      rules.alpha_div.output,
      rules.nmds.output,
      rules.nmds_hull.output,
      rules.correlation_clustering.output
   output:
      "{sample}.final"
   shell:
      "touch {output}"
#input: lambda wildcards: glob('sample-{samp}/*.fastq'.format(samp=wildcards.samp))