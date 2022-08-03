checkpoint biom_pairwise: # Create pairwise biom files for all the variables so that the differential abundance tests would work
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory("data/diff/pairwise–{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}/"))
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}"
   shell:
      " mkdir -p {output} && "
      " Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o {output} "

rule testDA:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/diff/pairwise–{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.tsv",
      map="data/diff/pairwise–{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.txt"
   output:
      testfile="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diff–{annie}–{theTest}.tsv"
   log: 
      "data/logs/DAtest–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Test {wildcards.theTest} for {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_test.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -o {output.testfile} -t {wildcards.theTest} > {log} 2>&1"


def ids_testDA(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diff–{annie}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names
                         
rule testDA_done: 
   input:
      ids_testDA
   output:
      temporary(touch("tmp/testDA_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))



rule EffSizePowerTest:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/diff/pairwise–{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.tsv",
      map="data/diff/pairwise–{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.txt"
   output:
      testfile="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/EffSizePowerTest–{annie}–{theTest}.tsv"
   log: 
      "data/logs/EffSizePowerTest–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Effect size power testing: Test {wildcards.theTest} for {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/powerDA.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -o {output.testfile} -t {wildcards.theTest} > {log} 2>&1"


def ids_EffSizePowerTest(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/EffSizePowerTest–{annie}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names
                         
rule EffSizePowerTest_done: 
   input:
      ids_EffSizePowerTest
   output:
      temporary(touch("tmp/EffSizePowerTest{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))













rule PowerDA:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/diff/pairwise–{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.tsv",
      map="data/diff/pairwise–{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}/{annie}.txt"
   output:
      testfile=temporary(directory("data/diff/Power_{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–{theTest}/"))
   log: 
      "data/logs/DAtest–Power–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Power analysis for the test {wildcards.theTest} using {wildcards.sample} – {wildcards.annie}."
   shell:
      "mkdir -p {output} && "
      "Rscript --vanilla workflow/scripts/testDA.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -t {wildcards.theTest} -g {output.testfile}/diffPower–{wildcards.annie}–{wildcards.theTest}.svg -o {output.testfile}/diffPower–{wildcards.annie}–{wildcards.theTest}.tsv > {log} 2>&1"


# Remove the annie and theTest from here on. Give it a try
def ids_PowertestDA(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise.get(**wildcards).output[0]    
    file_names = expand("data/diff/Power_{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–{theTest}/", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names

rule testDA_Power_plot: 
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      ids_PowertestDA
   output:
      plotfile="data/plots/DAscores–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      mytable ="data/diff/DAscores–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   log: 
      "data/logs/DAscores–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Plotting power analysis using {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_plot.R -i data/diff/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/ -g {output.plotfile} -o {output.mytable} > {log} 2>&1"

rule testDA_Power_test:
   conda:
      "../../workflow/envs/dunn.yaml"
   input:
      mytable ="data/diff/DAscores–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   output:
      dunnfile="data/plots/DAscores–Dunn–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.tsv"
   log: 
      "data/logs/DAscoresDunn–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Nonparametric testing (Dunn's multiple comparisons) power analysis using {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_dunn.R -i {input.mytable} -o {output.dunnfile} > {log} 2>&1"



rule testDAPower_done: 
   input:
      ids_PowertestDA
   output:
      temporary(touch("tmp/testDA_Power_{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))









rule diff:
   input:
      expand(         "tmp/testDA_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt", sample=config["mysample"],group=config["group"],theTest=config["DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"]),
      expand(   "tmp/testDA_Power_{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–done.txt", sample=config["mysample"],group=config["group"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"]),
      expand("tmp/EffSizePowerTest{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt", sample=config["mysample"],group=config["group"],theTest=config["DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"])

