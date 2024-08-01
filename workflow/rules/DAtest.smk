rule testDA:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.tsv",
      map="data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.txt"
   output:
      testfile=report("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diff–{annie}–{theTest}.tsv",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group} – {annie}",
      labels={
              "Description": "{theTest} – Differential abundance for samples",
              "Data type": "Text file"})
   log: 
      "data/logs/testDA–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Test {wildcards.theTest} for {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_test.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -o {output.testfile} -t {wildcards.theTest} > {log} 2>&1"





def ids_testDA(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise2.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diff–{annie}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names

rule testDA_done: 
   input:
      ids_testDA
   output:
      temporary(touch("data/tmp/testDA_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))














rule PowerDA:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.tsv",
      map="data/tmp/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.txt"
   output:
      testfile=report("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPower–{annie}–{theTest}.tsv",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group} – {annie}",
      labels={
              "Description": "{theTest} Spike-based Effect size",
              "Data type": "Text file"})
   log: 
      "data/logs/PowerDA–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Power analysis for the test {wildcards.theTest} using {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/powerDA.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -t {wildcards.theTest} -o {output.testfile} -l {log} > {log} 2>&1"



def ids_testDAPower(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise2.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPower–{annie}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names


rule testDAPower_done: 
   input:
      ids_testDAPower
   output:
      temporary(touch("data/tmp/testDA_Power_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))





rule cleanup_DAtest_pairwise:
   input:
      "data/tmp/testDA_Power_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt",
      "data/tmp/testDA_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"
   output:
      temporary(touch("data/tmp/done_{sample}–{group}–{theTest}–{minabund}{minread}{minsample}.txt"))
   message: "Done with testing effect size power testing for {wildcards.sample}–{wildcards.group}–{wildcards.theTest}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}"
   shell:
      "echo done"
















rule EffSizePowerTest:
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/tmp/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/pairwise/{annie}.tsv",
      map="data/tmp/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/pairwise/{annie}.txt"
   output:
      thefile=temporary("data/tmp/done–DATest–EffSizePowerTest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.txt")
   params:
      DAtest=config["DA_tests"]
   message: "Effect size power testing for {wildcards.sample} – {wildcards.annie}."
   shell:
      '''
      mkdir -p data/diff/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}/AUC_FDR_Power && echo {params.DAtest} | 
      tr " " "\n" | 
      parallel -j 1 "if [ ! -e data/diff/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}/AUC_FDR_Power/EffSizePowerTest-{{}}.tsv ] ; then Rscript --vanilla workflow/scripts/testDA.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -p /tmp/ -l data/logs/EffSizePowerTest–{wildcards.sample}–{wildcards.annie}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}–{{}}.txt -o data/diff/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}/AUC_FDR_Power/EffSizePowerTest-{{}}.tsv -t {{}} > data/logs/EffSizePowerTest–{wildcards.sample}–{wildcards.annie}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}–{{}}.txt  2>&1 && rm -rf /data/diff/pairwise–{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}–{{}}/ ;fi" &&
       touch {output.thefile}
      '''



rule testDA_Power_plot: 
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      thedonefile="data/tmp/done–DATest–EffSizePowerTest–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.txt"
   output:
      mytable ="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–AUC_FDR_Power.txt",
      myAUC=report("data/plots/AUC–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group} – {annie}",
      labels={
              "Description": "Area under the curve for differential abundance tests",
              "Data type": "svg file"}),
      myPower=report("data/plots/Power–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group} – {annie}",
      labels={
              "Description": "Power for differential abundance tests",
              "Data type": "svg file"}),
      myFDR=report("data/plots/FDR–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group} – {annie}",
      labels={
              "Description": "False discovery rate for differential abundance tests",
              "Data type": "svg file"}),
      myScore=report("data/plots/Score–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group} – {annie}",
      labels={
              "Description": "Composite test scores for differential abundance",
              "Data type": "svg file"})
   log: 
      "data/logs/testDA_Power_plot{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Plotting power analysis using {wildcards.sample} – {wildcards.annie}."
   shell:
      "mkdir -p ./data/plots | Rscript --vanilla workflow/scripts/DA_plot.R -i data/diff/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}/AUC_FDR_Power -o {output.mytable} -a {output.myAUC} -p {output.myPower} -f {output.myFDR} -s {output.myScore} > {log} 2>&1 "

rule testDA_Power_test: 
   conda:
      "../../workflow/envs/dunn.yaml"
   input:
      mytable ="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–AUC_FDR_Power.txt"
   output:
      dunnfile=report("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–Dunn.tsv",
      caption="../report/DAtest_subject_test.rst",
      category="Differential abundance – {sample} – {group} – {annie}",
      labels={
              "Description": "Dunn's nonparametric test for the scores",
              "Data type": "Text file"})
   log: 
      "data/logs/DAscoresDunn–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{annie}.txt"
   message: "Differential abundance: Dunn's multiple comparisons) power analysis using {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_dunn.R -i {input.mytable} -o {output.dunnfile} > {log} 2>&1 "

def ids_EffSizePowerTest(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise1.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–Dunn.tsv", sample=wildcards.sample, group=wildcards.group, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names


rule EffSizePowerTest_done: 
   input:
      ids_EffSizePowerTest
   output:
      temporary(touch("data/tmp/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–TestDA_power_done.txt"))
   message: "Done with power testing, graphing and Dunn's test {wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}"
   shell:
      "echo done"

rule diff:
   input:
      expand("data/tmp/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–TestDA_power_done.txt", sample=config["mysample"],group=config["group"],theTest=config["DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"]),
      expand("data/tmp/done_{sample}–{group}–{theTest}–{minabund}{minread}{minsample}.txt", sample=config["mysample"],group=config["group"],theTest=config["DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"])
