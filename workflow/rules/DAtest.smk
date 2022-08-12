def get_mem_mb(wildcards, attempt):
    return attempt * 1000


checkpoint biom_pairwise1: # Create pairwise biom files for all the variables so that the differential abundance tests would work
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/pairwise"))
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}"
   shell:
      " mkdir -p {output} && "
      " Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o {output} "

checkpoint biom_pairwise2: # Create pairwise biom files for all the variables so that the differential abundance tests would work
   conda:
      "../../workflow/envs/phyloseq_vegan_tidyverse.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      temporary(directory("data/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/"))
   message: "Filtering biom files of {wildcards.sample} as pairwise comparisons of the variables of {wildcards.group}"
   shell:
      " mkdir -p {output} && "
      " Rscript --vanilla workflow/scripts/pairwise_biom.R -i {input} -m data/{wildcards.sample}.txt -g {wildcards.group} -o {output} "


rule testDA:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.tsv",
      map="data/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.txt"
   output:
      testfile="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diff–{annie}–{theTest}.tsv"
   resources:
      mem_mb=get_mem_mb
   log: 
      "data/logs/DAtest–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
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
      temporary(touch("tmp/testDA_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))






rule PowerDA:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.tsv",
      map="data/diff/pairwise–{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{theTest}/{annie}.txt"
   output:
      testfile="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPower–{annie}–{theTest}.tsv"
   resources:
      mem_mb=get_mem_mb
   log: 
      "data/logs/DAtest–Power–{sample}–{group}–{annie}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}.txt"
   message: "Differential abundance: Power analysis for the test {wildcards.theTest} using {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/testDA.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -t {wildcards.theTest} -g {output}.svg -o {output.testfile} > {log} 2>&1"


def ids_testDAPower(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise2.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/diffPower–{annie}–{theTest}.tsv", sample=wildcards.sample, group=wildcards.group, theTest=wildcards.theTest, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names


rule testDAPower_done: 
   input:
      ids_testDAPower
   output:
      temporary(touch("tmp/testDA_Power_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"))





rule cleanup_DAtest_pairwise:
   input:
      "tmp/testDA_Power_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt",
      "tmp/testDA_{sample}–{group}–{theTest}–minAbd{minabund}minR{minread}minS{minsample}–done.txt"
   output:
      touch("tmp/done_{sample}–{group}–{theTest}–{minabund}{minread}{minsample}.txt")
   shell:
      "rm -rf data/diff/pairwise–{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}–{wildcards.theTest}/"





rule EffSizePowerTest:
   version: "1.0"
   conda:
      "../../workflow/envs/DAtest.yaml"
   input:
      tsv="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/pairwise/{annie}.tsv",
      map="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/pairwise/{annie}.txt"
   output:
      testfile=temporary(directory("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}/AUC_FDR_Power"))
   log: 
      "data/logs/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{annie}–AUC_FDR_Power.txt"
   resources:
      mem_mb=get_mem_mb
   params:
      DAtest=config["DA_tests"]
   message: "Effect size power testing for {wildcards.sample} – {wildcards.annie}."
   shell:
      '''
      mkdir -p {output} && echo {params.DAtest} | 
      tr " " "\n" | 
      parallel "echo Starting AUC, FDR and Power calculation of test {{}} | Rscript --vanilla workflow/scripts/testDA.R -i {input.tsv} -m {input.map} -c {wildcards.group} -s {wildcards.minsample} -r {wildcards.minread} -a {wildcards.minabund} -o {output.testfile}/EffSizePowerTest-{{}}.tsv -t {{}} >> {log} 2>&1"
      '''


rule testDA_Power_plot: 
   conda:
      "../../workflow/envs/ggpubr.yaml"
   input:
      "data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}/AUC_FDR_Power"
   output:
      plotfile="data/plots/AUC_FDR_Power–{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}.svg",
      mytable ="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–AUC_FDR_Power.txt"
   log: 
      "data/logs/{sample}–{group}–{annie}–minAbd{minabund}minR{minread}minS{minsample}–Plot.txt"
   message: "Differential abundance: Plotting power analysis using {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_plot.R -i {input} -g {output.plotfile} -o {output.mytable} > {log} 2>&1 && "
      "rm -rf data/diff/{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/{wildcards.annie}/"

rule testDA_Power_test: 
   conda:
      "../../workflow/envs/dunn.yaml"
   input:
      mytable ="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–AUC_FDR_Power.txt"
   output:
      dunnfile="data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–Dunn.tsv" 
   log: 
      "data/logs/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–{annie}–DAscoresDunn.txt"
   message: "Differential abundance: Nonparametric testing (Dunn's multiple comparisons) power analysis using {wildcards.sample} – {wildcards.annie}."
   shell:
      "Rscript --vanilla workflow/scripts/DA_dunn.R -i {input.mytable} -o {output.dunnfile} > {log} 2>&1"

def ids_EffSizePowerTest(wildcards):
    biom_pairwise_checkpoint_output = checkpoints.biom_pairwise1.get(**wildcards).output[0]    
    file_names = expand("data/diff/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}/{annie}–Dunn.tsv", sample=wildcards.sample, group=wildcards.group, minabund=wildcards.minabund, minread=wildcards.minread, minsample=wildcards.minsample, annie = glob_wildcards(os.path.join(biom_pairwise_checkpoint_output, "{annie}.tsv")).annie)
    return file_names


rule EffSizePowerTest_done: 
   input:
      ids_EffSizePowerTest
   output:
      touch("tmp/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–TestDA_power_done.txt")
   shell:
      "rm -rf data/diff/pairwise–{wildcards.sample}–{wildcards.group}–minAbd{wildcards.minabund}minR{wildcards.minread}minS{wildcards.minsample}/"

rule diff:
   input:
      expand("tmp/done_{sample}–{group}–{theTest}–{minabund}{minread}{minsample}.txt", sample=config["mysample"],group=config["group"],theTest=config["DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"]),
      expand("tmp/{sample}–{group}–minAbd{minabund}minR{minread}minS{minsample}–TestDA_power_done.txt", sample=config["mysample"],group=config["group"],theTest=config["DA_tests"],minabund=config["DA_minabund"],minread=config["DA_minread"],minsample=config["DA_minsample"])
      