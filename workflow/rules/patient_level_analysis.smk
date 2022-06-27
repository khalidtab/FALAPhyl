def get_mem_mb(wildcards, attempt):
    return attempt * 1000

rule patient_level_alpha_comparison:
   version: "1.0"
   conda: "../../workflow/envs/ggpubr.yaml"
   input:
      rules.alpha_div_calc.output
   output: 
      plots=report(directory("data/plots/patientlvl_alphaDiv_{sample}–{group}–{alpha}/"),
      caption="../report/ptlvl_alpha.rst",
      category="Patient-level analysis",
      patterns=["{name}.svg"],
      subcategory="Alpha – {alpha}",
      labels={
              "Data type": "Plot",
              "Comparison": "{name}",
              "Method": "{alpha}",
              "Grouping category": "{group}"})
   params:
      subjectID=expand("subjectID",subjectID=config["subjectID"])
   log:
      "data/logs/ptlevel_alpha–{sample}–{group}–{alpha}.txt"
   message: "{wildcards.alpha} Alpha diversity Patient-level calculations for {wildcards.sample}"
   shell:
      " mkdir {output} && "
      "Rscript --vanilla ./workflow/scripts/alpha_div_subjectlevel_nonparam.R -i {input} -g {wildcards.group} -m data/{wildcards.sample}.txt -p {params.subjectID} -o {output} > {log} 2>&1 "

      
      
rule patient_level_betapart: 
   version: "1.0"
   conda: "../../workflow/envs/ggpubr.yaml"
   input:
      repl=rules.betapart_matrix.output.repl,
      norepl=rules.betapart_matrix.output.norepl,
      jac=rules.betapart_matrix.output.bray
   resources:
      mem_mb=get_mem_mb
   output: 
      plots=report(directory("data/plots/patientlvl_betapart_{sample}–{group}–{distance}/"),
      caption="../report/ptlvl_beta.rst",
      category="Patient-level analysis",
      patterns=["{name}.svg"],
      subcategory="Beta – {distance} breakdown",
      labels={
              "Data type": "Plot",
              "Comparison": "{name}",
              "Method": "{distance}",
              "Grouping category": "{group}"})
   params:
      subjectID=expand("subjectID",subjectID=config["subjectID"])
   log:
      "data/logs/ptlevel_beta–{sample}–{distance}–{group}.txt"
   message: "Betapart {wildcards.distance} Patient-level calculations for {wildcards.sample}"
   shell:
      " mkdir {output} && "
      " Rscript --vanilla ./workflow/scripts/betapart_subjectlevel.R -i {input.jac} -r {input.repl} -n {input.norepl} -m data/{wildcards.sample}.txt -c {wildcards.group} -p {params.subjectID} -o {output} -d {wildcards.distance} > {log} 2>&1 "



rule subject_beta:
   input:
      expand("data/plots/patientlvl_betapart_{sample}–{group}–{distance}/",  sample=config["mysample"],group=config["group"],distance=["bray","jaccard"])

rule subject_alpha:
   input:
      expand("data/plots/patientlvl_alphaDiv_{sample}–{group}–{alpha}/",  sample=config["mysample"],group=config["group"],alpha=config["alpha"])