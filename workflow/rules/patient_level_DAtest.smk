rule testDA_pair with:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.tsv",
      map="data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.txt"
   output:
      testfile=report("data/diff_subject/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPaired–{annie}–{theTest}.tsv",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – Patient-level analysis – {sample} – {group} – {annie}",
      labels={
              "Description": "Test: {theTest}, Differential abundance for samples",
              "Data type": "Text file"})
   params:
      subject=config["subjectID"][0]
   log: 
      "data/logs/testDAPair–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Test {wildcards.theTest} for {wildcards.sample} – {wildcards.annie}, with within subject constraints."
   shell:
      "Rscript --vanilla workflow/scripts/DA_test_subject.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -o {output.testfile} -t {wildcards.theTest} -p {params.subject} > {log} 2>&1"



def ids_testDApair(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise2.get(**wildcards).output[0]    
    file_names = expand("data/diff_subject/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPaired–{annie}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names

rule testDAPair_done: 
   input:
      ids_testDApair
   output:
      temporary(touch("data/tmp/testDAPair_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))














rule PowerDAPair:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.tsv",
      map="data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.txt"
   output:
      testfile=report("data/diff_subject/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPowerPair–{annie}–{theTest}.tsv",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – Patient-level analysis – {sample} – {group} – {annie}",
      labels={
              "Description": "Test: {theTest}, Spike-based Effect size",
              "Data type": "Text file"})
   params:
      subject=config["subjectID"][0]
   log: 
      "data/logs/PowerDAPair–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Power analysis for the test {wildcards.theTest} using {wildcards.sample} – {wildcards.annie}, with within subject constraints."
   shell:
      "Rscript --vanilla workflow/scripts/powerDA_subject.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -t {wildcards.theTest} -o {output.testfile} -p {params.subject} -e /data/tmp > {log} 2>&1"

def ids_testDAPowerPair(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise2.get(**wildcards).output[0]    
    file_names = expand("data/diff_subject/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPowerPair–{annie}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names


rule testDAPowerPair_done: 
   input:
      ids_testDAPowerPair
   output:
      temporary(touch("data/tmp/testDAPair_Power_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))





rule cleanup_DAtestPair_pairwise:
   input:
      "data/tmp/testDAPair_Power_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt",
      "data/tmp/testDAPair_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"
   output:
      temporary(touch("data/tmp/donePair_{sample}–{group}–{theTest}–{minabund}{minread}{minsample}.txt"))
   message: "Done with testing effect size power testing, with within subject constraints, for {wildcards.sample}–{wildcards.group}–{wildcards.theTest}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}"

















rule EffSizePowerTestPair:
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/tmp/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/pairwise/{annie}.tsv",
      map="data/tmp/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/pairwise/{annie}.txt"
   output:
      thefile=temporary("data/tmp/done–DATestPair–EffSizePowerTest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.txt"),
   params:
      DAtest=config["paired_DA_tests"],
      subject=config["subjectID"][0]
   message: "Effect size power testing for {wildcards.sample} – {wildcards.annie}, with within subject constraints."
   shell:
      '''
      mkdir -p data/diff_subject/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}–Paired/AUC_FDR_Power && echo {params.DAtest} | 
      tr " " "\n" | 
      parallel -j 1 "if [ ! -e data/diff_subject/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}–Paired/AUC_FDR_Power/EffSizePowerTest-{{}}.tsv ] ; then Rscript --vanilla workflow/scripts/testDA_subject.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -o data/diff_subject/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}–Paired/AUC_FDR_Power/EffSizePowerTest-{{}}.tsv -t {{}} > data/logs/EffSizePowerTest–{wildcards.sample}–{wildcards.annie}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}–{{}}.txt -p {params.subject} 2>&1 ;fi" &&
       touch {output.thefile}
      '''


rule testDA_Power_plotPaired: 
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      thedonefile="data/tmp/done–DATestPair–EffSizePowerTest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.txt"
   output:
      mytable ="data/diff_subject/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–AUC_FDR_Power–Paired.txt",
      myAUC=report("data/plots/Paired–AUC–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – Patient-level analysis – {sample} – {group} – {annie}",
      labels={
              "Description": "Area under the curve for differential abundance tests",
              "Data type": "svg file"}),
      myPower=report("data/plots/Paired–Power–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – Patient-level analysis – {sample} – {group} – {annie}",
      labels={
              "Description": "Power for differential abundance tests",
              "Data type": "svg file"}),
      myFDR=report("data/plots/Paired–FDR–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – Patient-level analysis – {sample} – {group} – {annie}",
      labels={
              "Description": "False discovery rate for differential abundance tests",
              "Data type": "svg file"}),
      myScore=report("data/plots/Paired–Score–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – Patient-level analysis – {sample} – {group} – {annie}",
      labels={
              "Description": "Composite test scores for differential abundance",
              "Data type": "svg file"})
   log: "data/logs/testDAPair_Power_plot{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Plotting power analysis using {wildcards.sample} – {wildcards.annie}, with within subject constraints."
   shell:
      "mkdir -p ./data/plots | Rscript --vanilla workflow/scripts/DA_plot.R -i data/diff_subject/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}–Paired/AUC_FDR_Power -o {output.mytable} -a {output.myAUC} -p {output.myPower} -f {output.myFDR} -s {output.myScore} > {log} 2>&1 "


use rule testDA_Power_test as testDAPair_Power_test with: 
   input:
      mytable ="data/diff_subject/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–AUC_FDR_Power–Paired.txt"
   output:
      dunnfile=report("data/diff_subject/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–Dunn–Paired.tsv",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – Patient-level analysis – {sample} – {group} – {annie}",
      labels={
              "Description": "Dunn's nonparametric test for the scores",
              "Data type": "Text file"}) 
   log: 
      "data/logs/DAscoresDunnPaired–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{annie}.txt"
   message: "Differential abundance: Dunn's nonparametric test for  power analysis using {wildcards.sample} – {wildcards.annie}, with within subject constraints."





def ids_EffSizePowerTestPair(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise1.get(**wildcards).output[0]    
    file_names = expand("data/diff_subject/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–Dunn–Paired.tsv", sample=wildcards.sample, group=wildcards.group, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names


rule EffSizePowerTestPair_done: 
   input:
      ids_EffSizePowerTestPair
   output:
      temporary(touch("data/tmp/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–TestDAPair_power_done.txt"))
   message: "Done with power testing, graphing and Dunn's test, with within subject constraints. {wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}"
   shell:
      "echo done"


rule subject_diff:
   input:
      expand("data/tmp/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–TestDAPair_power_done.txt", sample=config["mysample"],group=config["group"],theTest=config["paired_DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"]),
      expand("data/tmp/donePair_{sample}–{group}–{theTest}–{minabund}{minread}{minsample}.txt", sample=config["mysample"],group=config["group"],theTest=config["paired_DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"])
