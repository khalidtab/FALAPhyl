def get_mem_mb(wildcards, attempt):
    return attempt * 1000

rule patient_level_alpha_comparison:
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
      subjectID=expand("{subjectID}",subjectID=config["subjectID"])
   log:
      "data/logs/ptlevel_alpha–{sample}–{group}–{alpha}.txt"
   message: "{wildcards.alpha} Alpha diversity Patient-level calculations for {wildcards.sample}"
   shell:
      " mkdir {output} && "
      "Rscript --vanilla ./workflow/scripts/alpha_div_subjectlevel_nonparam.R -i {input} -g {wildcards.group} -m data/{wildcards.sample}.txt -p {params.subjectID} -o {output} > {log} 2>&1 "

      
      
rule patient_level_betapart: 
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
      subjectID=expand("{subjectID}",subjectID=config["subjectID"])
   log:
      "data/logs/ptlevel_beta–{sample}–{distance}–{group}.txt"
   message: "Betapart {wildcards.distance} Patient-level calculations for {wildcards.sample}"
   shell:
      " mkdir {output} && "
      " Rscript --vanilla ./workflow/scripts/betapart_subjectlevel.R -i {input.jac} -r {input.repl} -n {input.norepl} -m data/{wildcards.sample}.txt -c {wildcards.group} -p {params.subjectID} -o {output} -d {wildcards.distance} > {log} 2>&1 "

rule friedman: # Plots the beta diversity distances using the Non-Metric Dimensional Scaling (NMDS) algorithm
   conda: "../../workflow/envs/dunn.yaml"
   input:
      rules.beta_div.output.tsv
   output: 
      report("data/distance/PerGroupDistCompare/{sample}/friedman–{dist}–{group}.txt",
      category="Beta diversity",
      subcategory="{dist}",
      labels={
              "Description": "Friedman's nonparametric test and Dunn's multiple comparisons of the distances between the groups",
              "Data type": "Text file",
              "Distance type": "{dist}",
              "Grouping category": "{group}"})
   log:
      "data/logs/Friendman_Dunn_{sample}–{dist}–{group}.log"
   message: "Beta diversity - {wildcards.dist}: Calculating distances between groups based on Friedman's nonparametric test, then performing Dunn's multiple comaparison for variable {wildcards.group} in {wildcards.sample}"
   params:
      subjectID=expand("{subjectID}",subjectID=config["subjectID"])
   shell:
      "Rscript --vanilla ./workflow/scripts/friedman.R -i {input} -o {output} -m data/{wildcards.sample}.txt -c {wildcards.group} -s {params.subjectID} > {log} 2>&1 "

rule paired_beta:
   input:
      expand("data/plots/patientlvl_betapart_{sample}–{group}–{distance}/",   sample=config["mysample"],group=config["group"],distance=["bray","jaccard"]),
      expand("data/distance/PerGroupDistCompare/{sample}/friedman–{dist}–{group}.txt",sample=config["mysample"],group=config["group"],dist=config["distances"])

rule paired_alpha:
   input:
      expand("data/plots/patientlvl_alphaDiv_{sample}–{group}–{alpha}/",  sample=config["mysample"],group=config["group"],alpha=config["alpha"])
