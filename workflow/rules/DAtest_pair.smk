rule testDApair with:
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:  # Input is from makeBiomForDAtest
      tsv="data/tmp/diff/DAtest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{pairwiseCats}.tsv",
      map="data/tmp/diff/DAtest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{pairwiseCats}.txt"
   output:
      testfile=report("data/diffpaired/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPaired–{pairwiseCats}–{theTest}.tsv",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – Patient-level analysis – {sample} – {group} – {pairwiseCats}",
      labels={
              "Description": "{theTest} – Differential abundance for paired samples",
              "Data type": "Text file"})
   params:
      subject=config["subjectID"][0]
   log: 
      "data/logs/testDAPair–{sample}–{group}–{pairwiseCats}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance Paired: Test {wildcards.theTest} for {wildcards.sample} – {wildcards.pairwiseCats}"
   shell:
      '''
      Rscript --vanilla workflow/scripts/DApair_test.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -o {output.testfile} -t {wildcards.theTest} -u {params.subject} > {log} 2>&1
      '''


def ids_testDApair(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.makeBiomForDAtest.get(**wildcards).output[0]    
    file_names = expand("data/diffpaired/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPaired–{pairwiseCats}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, pairwiseCats = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{pairwiseCats}.tsv")).pairwiseCats)
    return file_names

rule testDAPair_done: 
   input:
      ids_testDApair
   output:
      temporary(touch("data/tmp/testDAPair_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))


rule EffSizePowerTestPaired:
   conda:
      "../../workflow/envs/DAtest.yaml"
   input: # Input from makeBiomForEffectSize
      tsv="data/tmp/diffPaired/EffectSize–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{pairwiseCats}.tsv",
      map="data/tmp/diffPaired/EffectSize–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{pairwiseCats}.txt"
   output:
      thefile="data/diffpaired/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/AUC_FDR_Power/effectSize–{effectSize}–{theTest}/{pairwiseCats}.tsv"
   params:
      subjectID=config["subjectID"][0]
   log:
      "data/logs/EffectSizePair–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–effectSize–{effectSize}–{theTest}–{pairwiseCats}.log"
   message: "Effect size power testing for {wildcards.sample} – {wildcards.pairwiseCats} – test {wildcards.theTest} – Effect size {wildcards.effectSize}"
   shell:
      '''
      Rscript --vanilla workflow/scripts/DApair_effectsize.R -i {input.tsv} -t {wildcards.theTest} -e {wildcards.effectSize} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -u {params.subjectID} -p data/tmp/ -o {output} > {log} 2>&1 
      '''


def ids_EffSizePairedPowerTest(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.makeBiomForEffectSizePaired.get(**wildcards).output[0]    
    file_names = expand("data/diffpaired/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/AUC_FDR_Power/effectSize–{effectSize}–{theTest}/{pairwiseCats}.tsv", sample=wildcards.sample, group=wildcards.group, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample,theTest=config["paired_DA_tests"], effectSize=config["DA_effectSize"], pairwiseCats = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{pairwiseCats}.tsv")).pairwiseCats)
    return file_names


rule DA_EffectSizePaired_plot: 
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      ids_EffSizePairedPowerTest
   output: 
      mytable ="data/diffpaired/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/AUC_FDR_Power.txt",
      myAUC=report("data/plots/AUCpaired–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance paired – {sample} – {group}",
      labels={
              "Description": "Area under the curve for differential abundance paired tests",
              "Data type": "svg file"}),
      myPower=report("data/plots/Powerpaired–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance paired – {sample} – {group}",
      labels={
              "Description": "Power for differential abundance paired tests",
              "Data type": "svg file"}),
      myFDR=report("data/plots/FDRpaired–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance paired – {sample} – {group}",
      labels={
              "Description": "False discovery rate for differential abundance paired tests",
              "Data type": "svg file"}),
      myScore=report("data/plots/Scorepaired–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance paired – {sample} – {group}",
      labels={
              "Description": "Composite test scores for differential abundance paired tests",
              "Data type": "svg file"})
   log: 
      "data/logs/EffectSizePaired_plot–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance Paired: Plotting power analysis using {wildcards.sample}"
   shell:
      "mkdir -p ./data/plots | Rscript --vanilla workflow/scripts/DA_plot.R -i data/diffpaired/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/AUC_FDR_Power -o {output.mytable} -a {output.myAUC} -p {output.myPower} -f {output.myFDR} -s {output.myScore} > {log} 2>&1 "




rule EffectSizePaired_StatTest: 
   conda:
      "../../workflow/envs/dunn.yaml"
   input:
      mytable ="data/diffpaired/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/AUC_FDR_Power.txt"
   output:
      dunnfile=report("data/diffpaired/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/nonparam_test_for_scores.txt",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance paired – {sample} – {group}",
      labels={
              "Description": "Dunn's nonparametric test for the scores",
              "Data type": "Text file"})
   log: 
      "data/logs/DAscoresDunn–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance Paired: Dunn's multiple comparisons) power analysis using {wildcards.sample}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_dunn.R -i {input.mytable} -o {output.dunnfile} > {log} 2>&1 "





rule paired_diff:
   input: #first for effect size and second is testDA
      expand("data/diffpaired/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/nonparam_test_for_scores.txt", sample=config["mysample"],group=config["group"],theTest=config["paired_DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"],effectSize=config["DA_effectSize"]),
      expand("data/tmp/testDAPair_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt", sample=config["mysample"],group=config["group"],theTest=config["paired_DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"])

