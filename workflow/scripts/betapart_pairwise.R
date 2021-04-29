#!/usr/local/bin/Rscript --vanilla

set.seed(1234)

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(require(tidyverse)))
suppressWarnings(suppressMessages(library(optparse)))	
suppressWarnings(suppressMessages(library(ggpubr)))	


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Distance matrix", metavar="Bray-Curtis Distance matrix that is going to be deconstructed to pairwise distances"),
  make_option(c("-r", "--repl"), type="character", default=NULL, help="Bray Curtis (replacement) Distance matrix", metavar="Repl Bray-Curtis Distance matrix that is going to be deconstructed to pairwise distances"),
  make_option(c("-n", "--norepl"), type="character", default=NULL, help="Bray Curtis (no replacement) Distance matrix", metavar="Repl Bray-Curtis Distance matrix that is going to be deconstructed to pairwise distances"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-c", "--category"), type="character", default=NULL, help="Category", metavar="Column in the mapping file that is the category to be compared"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder for all pairwise files", metavar="Output folder name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load variables
cond = opt$category # cond = "HMPbodysubsite"

# Load Bray Curtis
bray = opt$input # bray = "~/Desktop/bioinfo_snakemake/data/distance/beta_div/v35_oral_gs+braycurtis.tsv"
bray = suppressWarnings(suppressMessages(read_tsv(bray,progress = FALSE))) %>% as.matrix(.)
rownames(bray) = bray[,1]
bray = bray[, colnames(bray) != "X1"]

bray[upper.tri(bray)] = NA
bray = reshape2::melt(bray, na.rm = TRUE)
colnames(bray) = c("v1","v2","bray")

# Load Replacement Bray Curtis
repl = opt$input # repl = "~/Desktop/bioinfo_snakemake/data/distance/beta_div/v35_oral_gs+BrayRepl.tsv"
repl = suppressWarnings(suppressMessages(read_tsv(repl,progress = FALSE))) %>% as.matrix(.)
rownames(repl) = repl[,1]
repl = repl[, colnames(repl) != "X1"]

repl[upper.tri(repl)] = NA
repl = reshape2::melt(repl, na.rm = TRUE)
colnames(repl) = c("v1","v2","repl")

# Load No Replacement Bray Curtis
noRepl = opt$input # noRepl = "~/Desktop/bioinfo_snakemake/data/distance/beta_div/v35_oral_gs+BrayNoRepl.tsv"
noRepl = suppressWarnings(suppressMessages(read_tsv(noRepl,progress = FALSE))) %>% as.matrix(.)
rownames(noRepl) = noRepl[,1]
noRepl = noRepl[, colnames(noRepl) != "X1"]

noRepl[upper.tri(noRepl)] = NA
noRepl = reshape2::melt(noRepl, na.rm = TRUE)
colnames(noRepl) = c("v1","v2","noRepl")

# Load mapping file
map = opt$mapping # map = "~/Desktop/bioinfo_snakemake/data/map/v35_oral_gs.txt"
map = suppressMessages(read_tsv(map))

output = opt$output # output = "~/Desktop/bioinfo_snakemake/data/plots/"


# Create pairwise table with all distances
distance = merge(bray,repl, all = TRUE) %>% merge(.,noRepl, all = TRUE)

sampleIDs = map[,1]
myConds = which(colnames(map) == cond) %>% map[,.] %>% as.matrix(.)
myUniqueConds = myConds %>% unique(.)

reducedMap = cbind(sampleIDs,myConds)
colnames(reducedMap) = cbind("v1","cond1")

distance = merge(distance, reducedMap, by.x=c("v1"), by.y=c("v1")) %>% merge(., reducedMap, by.x=c("v2"), by.y=c("v1"))
colnames(distance) = c("v2","v1","bray","repl","noRepl","cond1","cond2")

myCombos = combn(myUniqueConds,2) %>% t(.) # Get all 2 way combinations of your variables
uniqueCond1 = myCombos[,1] %>% unique(.)

for (myCond1 in uniqueCond1) {
  
  table1 = subset(distance,cond1 == myCond1)
  colnames(table1) = c("v1","v2","bray","repl","noRepl","cond1","cond2")
  
  table2 = subset(distance,cond2 == myCond1)
  colnames(table2) = c("v1","v2","bray","repl","noRepl","cond2","cond1")
  
  myTable = rbind(table1,table2)
  
  rm(table1,table2)
  
  cond2s = subset(as.data.frame(myCombos),V1 == myCond1) %>% .[,2]
  
  for (myCond2 in cond2s) {
  table2 = subset(myTable,cond2 == myCond2)
  pvalue = wilcox.test(table2$repl,table2$noRepl, paired = TRUE) %>% .$p.value
  
  
  pd = data.frame(replacement = table2$repl, noreplacement = table2$noRepl)
  colnames(pd) = c("Turn-over","Difference in richness")
 
  myGraph = ggpaired(pd, cond1 = "Turn-over", cond2 = "Difference in richness",
           fill = "condition", palette = "jco", line.size=0.01, 
           title= (paste("Bray-Curtis breakdown between ",myCond1,"and",myCond2, "Pvalue <",pvalue))) 
   
  svg(paste0(output,myCond1,"_",myCond2,".svg"))
 print(myGraph)
  dev.off()
    }}





