def get_mem_mb(wildcards, attempt):
    return attempt * 1000

rule patient_level_alpha_comparison:
   version: "1.0"
   conda: "../../workflow/envs/ggpubr.yaml"
   input:
      rules.alpha_div_calc.output
   output: 
      plots=directory("data/plots/patientlevel_alphaDiv_{sample}/")
   params:
      group=expand("{group}",group=config["group"]),
      alpha=expand("{alpha}",alpha=config["alpha"]),
      subjectID=expand("subjectID",subjectID=config["subjectID"])
   message: "Alpha diversity Patient-level calculations for {wildcards.sample}"
   shell:
      "mkdir -p data/plots/patientlevel_alphaDiv_{wildcards.sample}/ &&"
      "echo 'for x in {params.group}; do for y in {params.alpha}; do "
      "Rscript --vanilla ./workflow/scripts/alpha_div_subjectlevel_nonparam.R -i data/alpha_div/calc_{wildcards.sample}+$y.txt -g $x -m data/map/{wildcards.sample}.txt -p {params.subjectID} -o data/plots/patientlevel_alphaDiv_{wildcards.sample}/ "
      "; done ; done' > tmp/patientlevel_alphaDiv_{wildcards.sample}.sh &&"
      "chmod +x tmp/patientlevel_alphaDiv_{wildcards.sample}.sh &&"
      "bash tmp/patientlevel_alphaDiv_{wildcards.sample}.sh"

rule patient_level_betapart_jaccard: # Plots the PhILR distances using the Non-Metric Dimensional Scaling (NMDS) algorithm, and adds a hull around the samples
   version: "1.0"
   conda: "../../workflow/envs/ggpubr.yaml"
   input:
      repl=rules.jaccard_betapart_matrix.output.repl,
      norepl=rules.jaccard_betapart_matrix.output.norepl,
      jac=rules.jaccard_betapart_matrix.output.jaccard
   output: 
      plots=directory("data/plots/patientlevel_betapart_jaccard_{sample}/")
   resources:
      mem_mb=get_mem_mb
   params:
      group=expand("{group}",group=config["group"]),
      subjectID=expand("subjectID",subjectID=config["subjectID"])
   message: "Betapart (Jaccard) Patient-level calculations for {wildcards.sample}"
   shell:
      "mkdir -p data/plots/patientlevel_betapart_jaccard_{wildcards.sample}/ &&"
      "echo 'for x in {params.group}; do "
      "Rscript --vanilla ./workflow/scripts/betapart_subjectlevel.R -i {input.jac} -r {input.repl} -n {input.norepl} -m data/map/{wildcards.sample}.txt -c $x -p {params.subjectID} -o data/plots/patientlevel_betapart_jaccard_{wildcards.sample}/ -d jaccard "
      "; done' > tmp/patientlevel_betapart_jaccard_{wildcards.sample}.sh &&"
      "chmod +x tmp/patientlevel_betapart_jaccard_{wildcards.sample}.sh &&"
      "bash tmp/patientlevel_betapart_jaccard_{wildcards.sample}.sh"
      
      
      
rule patient_level_betapart_bray: # Plots the PhILR distances using the Non-Metric Dimensional Scaling (NMDS) algorithm, and adds a hull around the samples
   version: "1.0"
   conda: "../../workflow/envs/ggpubr.yaml"
   input:
      repl=rules.bray_betapart_matrix.output.repl,
      norepl=rules.bray_betapart_matrix.output.norepl,
      jac=rules.bray_betapart_matrix.output.bray
   resources:
      mem_mb=get_mem_mb
   output: 
      plots=directory("data/plots/patientlevel_betapart_bray_{sample}/")
   params:
      group=expand("{group}",group=config["group"]),
      subjectID=expand("subjectID",subjectID=config["subjectID"])
   message: "Betapart (Bray-Curtis) Patient-level calculations for {wildcards.sample}"
   shell:
      "mkdir -p data/plots/patientlevel_betapart_bray_{wildcards.sample}/ &&"
      "echo 'for x in {params.group}; do "
      "Rscript --vanilla ./workflow/scripts/betapart_subjectlevel.R -i {input.jac} -r {input.repl} -n {input.norepl} -m data/map/{wildcards.sample}.txt -c $x -p {params.subjectID} -o data/plots/patientlevel_betapart_bray_{wildcards.sample}/ -d bray"
      "; done' > tmp/patientlevel_betapart_bray_{wildcards.sample}.sh &&"
      "chmod +x tmp/patientlevel_betapart_bray_{wildcards.sample}.sh &&"
      "bash tmp/patientlevel_betapart_bray_{wildcards.sample}.sh"
      
      
      

rule patient_level_analysis: # Last step from betapart. Cleans up temporary files.
    input:
        alpha=rules.patient_level_alpha_comparison.output,
        beta_jac=rules.patient_level_betapart_jaccard.output,
        beta_bray=rules.patient_level_betapart_bray.output
    output:
        touch("tmp/.done_patientlevelanalysis_{sample}.txt")
    shell:
        "echo Done with patient-level analysis."
