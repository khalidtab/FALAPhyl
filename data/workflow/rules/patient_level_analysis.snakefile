def get_mem_mb(wildcards, attempt):
    return attempt * 1000

rule patient_level_alpha_comparison:
   version: "1.0"
   conda: "../../workflow/envs/ggpubr.yaml"
   input:
      rules.alpha_div_calc.output
   output: 
      plots=directory("data/plots/patientlvl_alphaDiv_{sample}–{group}–{alpha}/")
   params:
      subjectID=expand("subjectID",subjectID=config["subjectID"])
   log:
      "data/logs/ptlevel_alpha–{sample}–{group}–{alpha}.txt"
   message: "{wildcards.alpha} Alpha diversity Patient-level calculations for {wildcards.sample}"
   shell:
      " mkdir {output} && "
      "Rscript --vanilla ./workflow/scripts/alpha_div_subjectlevel_nonparam.R -i {input} -g {wildcards.group} -m data/map/{wildcards.sample}.txt -p {params.subjectID} -o {output} > {log} 2>&1 "

      
      
rule patient_level_betapart: # Plots the PhILR distances using the Non-Metric Dimensional Scaling (NMDS) algorithm, and adds a hull around the samples
   version: "1.0"
   conda: "../../workflow/envs/ggpubr.yaml"
   input:
      repl=rules.betapart_matrix.output.repl,
      norepl=rules.betapart_matrix.output.norepl,
      jac=rules.betapart_matrix.output.bray
   resources:
      mem_mb=get_mem_mb
   output: 
      plots=directory("data/plots/patientlvl_betapart_{sample}–{group}–{distance}/")
   params:
      subjectID=expand("subjectID",subjectID=config["subjectID"])
   log:
      "data/logs/ptlevel_beta–{sample}–{distance}–{group}.txt"
   message: "Betapart {wildcards.distance} Patient-level calculations for {wildcards.sample}"
   shell:
      " mkdir {output} && "
      " Rscript --vanilla ./workflow/scripts/betapart_subjectlevel.R -i {input.jac} -r {input.repl} -n {input.norepl} -m data/map/{wildcards.sample}.txt -c {wildcards.group} -p {params.subjectID} -o {output} -d {wildcards.distance} > {log} 2>&1 "



rule subject_beta:
   input:
      expand("data/plots/patientlvl_betapart_{sample}–{group}–{distance}/",  sample=config["mysample"],group=config["group"],distance=["bray","jaccard"])

rule subject_alpha:
   input:
      expand("data/plots/patientlvl_alphaDiv_{sample}–{group}–{alpha}/",  sample=config["mysample"],group=config["group"],alpha=config["alpha"])