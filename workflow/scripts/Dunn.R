#!/usr/local/bin/Rscript --vanilla

suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(dunn.test)))
suppressWarnings(suppressMessages(library(magrittr)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Distance matrix", metavar="Dissimilarity matrix that is going to be deconstructed to pairwise dissimilarity"),
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

# Functions
load_and_fix_Dis = function(myInput){
  
  myDis = read.csv(myInput,comment="",sep = "\t",row.names = 1) 
  myDis[lower.tri(myDis)] = NA
  myDis[,1] = row.names(myDis)
  colnames(myDis)[1] = "X"
  myDis = suppressMessages(reshape2::melt(myDis,na.rm = TRUE))
  colnames(myDis) = c("SampleID","variable","value")
  myDis = myDis[!(myDis$SampleID == myDis$variable),] # Delete distances from self
  myDis = left_join(myDis,mymap,by="SampleID")
  colnames(myDis) = c("SampleID1","SampleID","value1","condition1")
  myDis = left_join(myDis,mymap,by="SampleID")
  colnames(myDis) = c("SampleID1","SampleID2","value1","condition1","condition2")
  myDis = cbind(myDis$condition1,myDis$condition2,myDis$value1)
  colnames(myDis) = c("condition1","condition2","value")
  myDis = as.data.frame(myDis)
  
  TorF = myDis$condition1 > myDis$condition2 # Make sure that category 1 and 2 are always in the right order so that we don't, for example, have an issue where the same 2 categories just switched placed in the two columns
  myDis = cbind(myDis,TorF)
  myDisTRUE = subset(myDis,TorF == "TRUE")
  myDisFALSE = subset(myDis,TorF == "FALSE")
  myDisFALSE = myDisFALSE %$% data_frame(condition2,condition1,value,TorF)
  colnames(myDisFALSE) = c("condition1","condition2","value","TorF")
  
  myDis = rbind(myDisFALSE,myDisTRUE)
  
  myDis = cbind(myDis,paste0(myDis$condition1,"+",myDis$condition2))
  colnames(myDis) = c("condition1","condition2","value","TorF","comparison")
  myDis = myDis %$% data.frame(comparison,value)
  
  return(myDis)
}
GiveMeDunn = function(myDis){
  myDunn = dunn.test(myDis$value %>% as.numeric(.),
                     myDis$comparison,
                     method="bh",kw=FALSE,table=FALSE,list=FALSE) %$% data.frame(comparisons,Z,P,P.adjusted)
  myDunn[c('Category 1', 'Category 2')] = str_split_fixed(myDunn$comparison, ' - ', 2)
  DistanceInCategory1BiggerThanCategory2 = (myDunn$Z > 0)
  myDunn = cbind(myDunn,DistanceInCategory1BiggerThanCategory2) 
  myDunn = myDunn %$% data.frame(`Category 1`,`Category 2`,Z,P,P.adjusted,DistanceInCategory1BiggerThanCategory2)
  return(myDunn)
}

# Load and format mapping file
mymap = read.table(opt$mapping,comment="@",header=TRUE)
colnames(mymap)[1] = "SampleID"
catNum = which(colnames(mymap) == opt$category)
mymap = cbind(mymap[,1],mymap[,catNum]) %>% as.data.frame(.)
colnames(mymap) = c("SampleID","condition")

myDis = load_and_fix_Dis(opt$input) # Load dissimilarity matrix, "melt" it, and make sure that category 1 and category 2 do not switch places
myStats = myDis%>% group_by(comparison)%>% summarise(Median=median(as.numeric(value)), Max=max(as.numeric(value)), Min=min(as.numeric(value)), IQR=IQR(as.numeric(value)))


myKruskall =  myDis %$% kruskal.test(value,comparison)

if (myKruskall$p.value < 0.05){
if (length(unique(myDis$comparison)) == 2){
sink(opt$output)
print("Only 2 categories are available. No Dunn's comparison can be done. Below is the Kruskal-Wallis results, and distance characteristics")  
myKruskall
myStats
sink()
}else{

sink(opt$output)
myKruskall
myStats
sink()
myDunn = GiveMeDunn(myDis)
myTable = left_join(myDunn,myStats, by=c("Category.1"="comparison"))
myTable = left_join(myTable,myStats, by=c("Category.2"="comparison"))
colnames(myTable) = c("Category.1","Category.2","Z","P","P.adjusted","DistanceInCategory1BiggerThanCategory2","Median.Category1","Max.Category1","Min.Category1","IQR.Category1","Median.Category2","Max.Category2","Min.Category2","IQR.Category2")
suppressWarnings(write.table(myTable,file = opt$output,append=TRUE,quote=FALSE,sep="\t",row.names = FALSE)) # Do Dunn test and format it as a table


}} else {
  if (length(unique(myDis$comparison)) == 2){
    sink(opt$output)
    print("Only 2 categories are available. No Dunn's comparison can be done. Below is the Kruskal-Wallis results")  
    myKruskall
    sink()
  }else{ sink(opt$output)
  print("Kruskal-Wallis Rank Sum Test p-value is larger than 0.05, as such, Dunn test below may not correct.")
  sink()
  suppressWarnings(write.table(GiveMeDunn(myDis),file = opt$output,append=TRUE,quote=FALSE,sep="\t",row.names = FALSE)) # Do Dunn test and format it as a table
}}
