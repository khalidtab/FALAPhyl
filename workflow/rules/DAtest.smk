rule testDA:
   conda:
      "../../workflow/envs/DAtest.yaml"
   input: # Input is from makeBiomForDAtest
      tsv="data/tmp/diff/DAtest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{pairwiseCats}.tsv",
      map="data/tmp/diff/DAtest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{pairwiseCats}.txt"
   output:
      testfile=report("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diff–{pairwiseCats}–{theTest}.tsv",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group} – {pairwiseCats}",
      labels={
              "Description": "{theTest} – Differential abundance for samples",
              "Data type": "Text file"})
   log: 
      "data/logs/testDA–{sample}–{group}–{pairwiseCats}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Test {wildcards.theTest} for {wildcards.sample} – {wildcards.pairwiseCats}."
   shell:
      '''
      Rscript --vanilla workflow/scripts/DA_test.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -o {output.testfile} -t {wildcards.theTest} > {log} 2>&1
      '''




def ids_testDA(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.makeBiomForDAtest.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diff–{pairwiseCats}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, pairwiseCats = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{pairwiseCats}.tsv")).pairwiseCats)
    return file_names

rule testDA_done: 
   input:
      ids_testDA
   output:
      temporary(touch("data/tmp/testDA_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))


rule EffSizePowerTest:
   conda:
      "../../workflow/envs/DAtest.yaml"
   input: # Input from makeBiomForEffectSize
      tsv="data/tmp/diff/EffectSize–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{pairwiseCats}.tsv",
      map="data/tmp/diff/EffectSize–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{pairwiseCats}.txt"
   output:
      thefile="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/AUC_FDR_Power/effectSize–{effectSize}–{theTest}/{pairwiseCats}.tsv"
   log:
      "data/logs/EffectSize–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–effectSize–{effectSize}–{theTest}–{pairwiseCats}.log"
   message: "Effect size power testing for {wildcards.sample} – {wildcards.pairwiseCats} – test {wildcards.theTest} – Effect size {wildcards.effectSize}"
   shell:
      '''
      Rscript --vanilla workflow/scripts/DA_effectsize.R -i {input.tsv} -t {wildcards.theTest} -e {wildcards.effectSize} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -p data/tmp/ -o {output} > {log} 2>&1 
      '''




def ids_EffSizePowerTest(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.makeBiomForEffectSize.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/AUC_FDR_Power/effectSize–{effectSize}–{theTest}/{pairwiseCats}.tsv", sample=wildcards.sample, group=wildcards.group, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample,theTest=config["DA_tests"], effectSize=config["DA_effectSize"], pairwiseCats = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{pairwiseCats}.tsv")).pairwiseCats)
    return file_names


rule DA_EffectSize_plot: 
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      ids_EffSizePowerTest
   output: 
      mytable ="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/AUC_FDR_Power.txt",
      myAUC=report("data/plots/AUC–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group}",
      labels={
              "Description": "Area under the curve for differential abundance tests",
              "Data type": "svg file"}),
      myPower=report("data/plots/Power–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group}",
      labels={
              "Description": "Power for differential abundance tests",
              "Data type": "svg file"}),
      myFDR=report("data/plots/FDR–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group}",
      labels={
              "Description": "False discovery rate for differential abundance tests",
              "Data type": "svg file"}),
      myScore=report("data/plots/Score–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group}",
      labels={
              "Description": "Composite test scores for differential abundance",
              "Data type": "svg file"})
   log: 
      "data/logs/EffectSize_plot–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Plotting power analysis using {wildcards.sample}"
   shell:
      "mkdir -p ./data/plots | Rscript --vanilla workflow/scripts/DA_plot.R -i data/diff/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/AUC_FDR_Power -o {output.mytable} -a {output.myAUC} -p {output.myPower} -f {output.myFDR} -s {output.myScore} > {log} 2>&1 "




rule EffectSize_StatTest: 
   conda:
      "../../workflow/envs/dunn.yaml"
   input:
      mytable ="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/AUC_FDR_Power.txt"
   output:
      dunnfile=report("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/nonparam_test_for_scores.txt",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group}",
      labels={
              "Description": "Dunn's nonparametric test for the scores",
              "Data type": "Text file"})
   log: 
      "data/logs/DAscoresDunn–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Dunn's multiple comparisons) power analysis using {wildcards.sample}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_dunn.R -i {input.mytable} -o {output.dunnfile} > {log} 2>&1 "





rule diff:
   input: #first for effect size and second is testDA
      expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/nonparam_test_for_scores.txt", sample=config["mysample"],group=config["group"],theTest=config["DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"],effectSize=config["DA_effectSize"]),
      expand("data/tmp/testDA_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt", sample=config["mysample"],group=config["group"],theTest=config["DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"])
