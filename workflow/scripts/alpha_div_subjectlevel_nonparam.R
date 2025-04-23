suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library("optparse")))
suppressPackageStartupMessages(suppressWarnings(library(reshape2)))
suppressWarnings(suppressMessages(library(effectsize)))

option_list = list(
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="The mapping file", metavar="Mapping file"),
  make_option(c("-i", "--input"), type="character", default=NULL, help="The alpha diversity file", metavar="Alpha diversity file"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="The groups column in your mapping file you are comparing", metavar="Groups"),
  make_option(c("-p", "--patientID"), type="character", default=NULL, help="Column that has the patient ID", metavar="Patient ID"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Path to output file of plots", metavar="Path to output folder")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

# Load files
map = opt$mapping
alpha = opt$input
group = opt$group
output = opt$output 
patientID = opt$patientID 

map = suppressMessages(read_tsv(map))
colnames(map)[1] = "SampleID"
alpha = suppressMessages(read_tsv(alpha))
alphaDivType = colnames(alpha)[2]
# Merge then keep only the pertinent columns
myNames = c("SampleID","alpha")
colnames(alpha) = myNames
map_first_column = colnames(map)[1]
merged = merge(alpha,map, by.x = "SampleID", by.y = colnames(map)[1]) %>% .[, c(myNames, group, patientID)]
colnames(merged)[3] = "group"

myCombos = unique(merged$group) %>% combn(.,2)

numOfComparisons = dim(myCombos)[2]

for (x in 1:numOfComparisons){
  category1 = myCombos[,x][1]
  category2 = myCombos[,x][2]
  category1Table = subset(as.data.frame(merged),group == category1)
  category2Table = subset(as.data.frame(merged),group == category2)
  currenTable = rbind(category1Table,category2Table)
  patientIDrepetitions = table(currenTable[paste(patientID)])
  if(length(unique(patientIDrepetitions)) > 1){
    print(paste("Categories to be compared are",category1,category2))
    print(patientIDrepetitions)
    print("Your design is not balanced. You have some samples missing that preclude analysis of your subject-level alpha diversity comparison.")
  }
  
  patientIDtwos = subset(as.data.frame(patientIDrepetitions),Freq == 2) %>% .[,1] %>% as.matrix(.)
  category1reps = rep(category1,length(patientIDtwos))
  category2reps = rep(category2,length(patientIDtwos))
  
  table1 = cbind(patientIDtwos,category1reps)
  colnames(table1) = c(patientID,"group")
  table1merged = merge(table1,merged,by.all = TRUE)
  
  table2 = cbind(patientIDtwos,category2reps)
  colnames(table2) = c(patientID,"group")
  table2merged = merge(table2,merged,by.all = TRUE)
  
  tablesMerged = merge(table1merged,table2merged,by = patientID)
  tablesMerged = cbind(tablesMerged$alpha.x,tablesMerged$alpha.y) %>% as.data.frame(.)
  colnames(tablesMerged) = c(category1,category2)
  
  if (dim(tablesMerged)[1] == 0) { # Meaning, if you merged two categories and you had zero rows because none of the samples matched, exit loop and continue
    break
  }
  
  #Otherwise, if there are matches, continue
  
  pvalue = wilcox.test(tablesMerged[,1],tablesMerged[,2], paired = TRUE) %>% .$p.value
  effectSize = rank_biserial(tablesMerged[,1],tablesMerged[,2], paired = TRUE) %>% interpret(., rules = "funder2019")
  numOfPatients = length(tablesMerged[,1])
  
  myGraph = ggpubr::ggpaired(tablesMerged, cond1 = category1, cond2 = category2, fill = "condition", palette = "jco", line.size=0.01, 
                             title= (paste("Patient-level comparison of",alphaDivType,"\nBetween",category1,"&",category2)), 
                             subtitle = paste("n =",numOfPatients,".Wilcoxon signed rank sum test, P-value <",formatC(pvalue, format = "e", digits = 5),
                             "\nRank biserial effect size= ",formatC(effectSize$r_rank_biserial, format = "g", digits = 2)))
  ggsave(filename=paste0(output,"/",category1,"_",category2,".svg"),plot=myGraph)
}

