#!/usr/local/bin/Rscript --vanilla

set.seed(1234)

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(require(tidyverse)))
suppressWarnings(suppressMessages(library(optparse)))	
suppressWarnings(suppressMessages(library(ggpubr)))	
suppressWarnings(suppressMessages(library(ggplot2)))	

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Distance matrix", metavar="Bray-Curtis or Jaccard Distance matrix that is going to be deconstructed to pairwise distances"),
  make_option(c("-d", "--distance"), type="character", default=NULL, help="Type of distance matrix. Either the string 'bray' or 'jaccard'.", metavar="Distance matrix"),
  make_option(c("-r", "--repl"), type="character", default=NULL, help="Bray Curtis (replacement) Distance matrix", metavar="Repl Bray-Curtis or Jaccard Distance matrix that is going to be deconstructed to pairwise distances"),
  make_option(c("-n", "--norepl"), type="character", default=NULL, help="Bray Curtis (no replacement) Distance matrix", metavar="Repl Bray-Curtis Distance or Jaccard matrix that is going to be deconstructed to pairwise distances"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-c", "--category"), type="character", default=NULL, help="Category", metavar="Column in the mapping file that is the category to be compared"),
  make_option(c("-p", "--patientID"), type="character", default=NULL, help="Column that has the patient ID", metavar="Patient ID"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder for all pairwise files", metavar="Output folder name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load variables
cond = opt$category
dist_type = opt$distance

patientID = opt$patientID

# Load distance
dist = opt$input 
dist = suppressWarnings(suppressMessages(read_tsv(dist,progress = FALSE))) %>% as.matrix(.)
rownames(dist) = dist[,1]
dist = dist[, colnames(dist) != "X1"]

dist[upper.tri(dist)] = NA
dist = reshape2::melt(dist, na.rm = TRUE)
colnames(dist) = c("v1","v2",dist_type)

# Load Replacement
repl = opt$repl 
repl = suppressWarnings(suppressMessages(read_tsv(repl,progress = FALSE))) %>% as.matrix(.)
rownames(repl) = repl[,1]
repl = repl[, colnames(repl) != "X1"]

repl[upper.tri(repl)] = NA
repl = reshape2::melt(repl, na.rm = TRUE)
colnames(repl) = c("v1","v2","repl")

# Load No Replacement
noRepl = opt$norepl 
noRepl = suppressWarnings(suppressMessages(read_tsv(noRepl,progress = FALSE))) %>% as.matrix(.)
rownames(noRepl) = noRepl[,1]
noRepl = noRepl[, colnames(noRepl) != "X1"]

noRepl[upper.tri(noRepl)] = NA
noRepl = reshape2::melt(noRepl, na.rm = TRUE)
colnames(noRepl) = c("v1","v2","noRepl")

# Load mapping file
map = opt$mapping 
map = suppressMessages(read_tsv(map))

output = opt$output 


# Create pairwise table with all distances
distance = merge(dist,repl, all = TRUE) %>% merge(.,noRepl, all = TRUE)

sampleIDs = map[,1]
patientIDcolumn = which(colnames(map) == patientID)
patientIDcolumn = map[,patientIDcolumn]
sampleIDtable = cbind(sampleIDs,patientIDcolumn)

myConds = which(colnames(map) == cond) %>% map[,.] %>% as.matrix(.)
myUniqueConds = myConds %>% unique(.)

reducedMap = cbind(sampleIDs,myConds)
colnames(reducedMap) = cbind("v1","cond1")

distance = merge(distance, reducedMap, by.x=c("v1"), by.y=c("v1")) %>% merge(., reducedMap, by.x=c("v2"), by.y=c("v1"))

colnames(sampleIDtable) = c("v1","v1subjectID")
distance = merge(distance,sampleIDtable, by.x=c("v1"), by.y=c("v1"))

colnames(sampleIDtable) = c("v2","v2subjectID")
distance = merge(distance,sampleIDtable, by.x=c("v2"), by.y=c("v2"))

# Subset to only those from the same patient
distance = distance[distance$v1subjectID==distance$v2subjectID, ]

colnames(distance) = c("v2","v1",dist_type,"repl","noRepl","cond1","cond2","v1subjectID","v2subjectID")

uniqueCond = c(distance$cond1, distance$cond2) %>% unique(.)
comboUniq = combn(uniqueCond,2,simplify = FALSE)

for (myConds in comboUniq) {
  myCond1 = myConds[1]
  myCond2 = myConds[2]
  
  table1 = subset(distance,cond1 == myCond1)
  table1 = data.frame(table1$v1,table1$v2,table1[,3],table1$repl,table1$noRepl,table1$cond1,table1$cond2) # Rearrange the columns
  colnames(table1) = c("v1","v2",dist_type,"repl","noRepl","cond1","cond2") # Fix the names of the columns
  
  table2 = subset(distance,cond2 == myCond1)
  colnames(table2) = c("v1","v2",dist_type,"repl","noRepl","cond2","cond1") # This makes table2 to have the same format as table1
  table2 = data.frame(table2$v1,table2$v2,table2[,3],table2$repl,table2$noRepl,table2$cond1,table2$cond2) # Rearrange the columns
  colnames(table2) = c("v1","v2",dist_type,"repl","noRepl","cond1","cond2") # Fix the names
  
  myTable = rbind(table1,table2)
  # By now myTable cond1 equals myCond1. Next, subset mytable so that cond2 equals myCond2
  myTable = subset(myTable,cond2 == myCond2)
  
  i = c(3,4,5)
  myTable[,i] = apply(myTable[,i],2, function(x) as.numeric(as.character(x)))
  rm(table1,table2)
  
  pvalue = wilcox.test(myTable$repl,myTable$noRepl, paired = TRUE) %>% .$p.value
  
  
  pd = data.frame(replacement = myTable$repl, noreplacement = myTable$noRepl)
  
  numOfParticipants = dim(pd)[1]
  
  if (dist_type == "bray"){
  colnames(pd) = c("balanced variation in abundance","abundance gradients")
  myGraph = ggpubr::ggpaired(pd, cond1 = "balanced variation in abundance", cond2 = "abundance gradients", fill = "condition", palette = "jco", line.size=0.01, subtitle = paste0("n =",numOfParticipants,". Wilcoxon signed rank sum test, P-value < ",formatC(pvalue, format = "e", digits = 5)), title= (paste(dist_type, "breakdown \n Between ",myCond1,"and",myCond2, "\n Pvalue <",formatC(pvalue, format = "e", digits = 2))))
  } else {
    colnames(pd) = c("Turn-over","Nestedness")
    myGraph = ggpubr::ggpaired(pd, cond1 = "Turn-over", cond2 = "Nestedness", fill = "condition", palette = "jco", line.size=0.01, title= (paste(dist_type, "breakdown \n Between ",myCond1,"and",myCond2)), subtitle = paste0("n =",numOfParticipants,"\n Wilcoxon signed rank sum test, P-value < ",formatC(pvalue, format = "e", digits = 5)))  
  } 
 
  ggsave(filename=paste0(output,dist_type,"ptlvl_",myCond1,"_",myCond2,".svg"),plot=myGraph) 
  
  }





