#!/usr/local/bin/Rscript --vanilla

suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(dunn.test)))
suppressWarnings(suppressMessages(library(magrittr)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Distance matrix", metavar="Dissimilarity matrix that is going to be deconstructed to pairwise dissimilarity"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-c", "--category"), type="character", default=NULL, help="Category", metavar="Column in the mapping file that is the category to be compared"),
  make_option(c("-s", "--subjectID"), type="character", default=NULL, help="Subject ID", metavar="column in mapping file that has the subject ID"),
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
  
  myDis = read.csv(myInput,comment="",sep = "\t",row.names = 1) %>% as.data.frame(.)
  myDis[lower.tri(myDis)] = NA
  myDis[,1] = row.names(myDis)
  colnames(myDis)[1] = "X"
  myDis = suppressMessages(reshape2::melt(myDis,na.rm = TRUE))
  colnames(myDis) = c("SampleID","SampleID2","value")
  myDis = myDis[!(myDis$SampleID == myDis$SampleID2),] # Delete distances from self
  myDis = left_join(myDis,mymap,by="SampleID") # Join by the first sample only
  colnames(myDis) = c("SampleID1","SampleID","value1","subjectID1","condition1")
  myDis = left_join(myDis,mymap,by="SampleID")
  colnames(myDis) = c("SampleID1","SampleID","value1","subjectID1","condition1","subjectID2","condition2")
  myDis = myDis[(myDis$subjectID1 == myDis$subjectID2),] # Delete distances from samples not from the same subject
  
  TorF = myDis$condition1 > myDis$condition2 # Make sure that category 1 and 2 are always in the right order so that we don't, for example, have an issue where the same 2 categories just switched placed in the two columns
  myDis = cbind(myDis,TorF)
  myDisTRUE = subset(myDis,TorF == "TRUE")
  colnames(myDisTRUE) = c("SampleID1","SampleID2","value1","subjectID1","condition1","subjectID2","condition2","TorF")
  myDisTRUE = myDisTRUE %$% data_frame(subjectID1,condition1,condition2,value1)
  colnames(myDisTRUE) = c("subjectID","condition1","condition2","value")
  
  myDisFALSE = subset(myDis,TorF == "FALSE")
  myDisFALSE = myDisFALSE %$% data_frame(subjectID1,condition2,condition1,value1)
  colnames(myDisFALSE) = c("subjectID","condition1","condition2","value")
  
  myDis = rbind(myDisFALSE,myDisTRUE)
  
  myDis = cbind(myDis,paste0(myDis$condition1,"+",myDis$condition2))
  colnames(myDis) = c("subjectID","condition1","condition2","value","comparison")
  myDis = myDis %$% data.frame(subjectID,comparison,value)
  
  
  
  return(myDis)
}
GiveMeDunn = function(myDis){
  mycolumns = colnames(myDis)

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
mymap = suppressMessages(read.csv(opt$mapping, skip=0, header=T, sep="\t"))
colnames(mymap)[1] = "SampleID"
catNum = which(colnames(mymap) == opt$category)
catSubject = which(colnames(mymap) == opt$subjectID)
mymap = cbind(mymap[,1],mymap[,catSubject],mymap[,catNum]) %>% as.data.frame(.)
colnames(mymap) = c("SampleID","subjectID","condition")

myDis = load_and_fix_Dis(opt$input) # Load dissimilarity matrix, "melt" it, and make sure that category 1 and category 2 do not switch places, and deletes distances from the same person, and any subjects that we don't have complete set of their sites

myStats = myDis %>% group_by(comparison)%>% summarise(Median=median(as.numeric(value)), Max=max(as.numeric(value)), Min=min(as.numeric(value)), IQR=IQR(as.numeric(value)))

myDis2 = pivot_wider(myDis,id_cols = "subjectID", names_from = "comparison", values_from = "value")

# At this point, you're done with processing the data. But, sometimes, when you don't have full datasets of the patients, you need to filter them out
myDis2 = na.omit(myDis2) %>% as.data.frame(.)
print(paste("The original number of samples is",dim(myDis2)[1],"and those that were retained after removal of samples without complete sets is",dim(myDis2)[1]))
rownames(myDis2) = myDis2[,1]
myDis2[,1] = NULL

myFriedman =  stats::friedman.test(as.matrix(myDis2))

if (myFriedman$p.value < 0.05){
if (length(unique(myDis$comparison)) == 2){

sink(opt$output)
print("Only 2 categories are available. No Dunn's comparison can be done. Below is the Friedman's nonparametric test results.")  
myFriedman
myStats
sink()
}else{

sink(opt$output)
myFriedman
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
    print("Only 2 categories are available. No Dunn's comparison can be done. Below is the Friedman's nonparametric test results")  
    myKruskall
    sink()
  }else{ sink(opt$output)
  print("Friedman's nonparametric Test p-value is larger than 0.05, as such, Dunn test below may not correct.")
  sink()
  suppressWarnings(write.table(GiveMeDunn(myDis),file = opt$output,append=TRUE,quote=FALSE,sep="\t",row.names = FALSE)) # Do Dunn test and format it as a table
}}
