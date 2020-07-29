library(optparse)
library(phyloseq)
suppressMessages(library(tidyverse))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input biom file", metavar="Input biom file"),
  make_option(c("-m", "--map"), type="character", default=NULL, help="The mapping file", metavar="Mapping file"),
  make_option(c("-c", "--categories"), type="character", default=NULL, help="The categories to filter by", metavar="Categories"),
  make_option(c("-t", "--threshold"), type="character", default=NULL, help="Threshold for core filtering", metavar="Threshold"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder for core files", metavar="output"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="the group from the mapping file to filter by", metavar="group")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

input = opt$input # input ="/Users/khaled/Desktop/bioinfo_snakemake/data/network/kfp5s2_CPd_GAPd_LAPd_corr/corr+CPd.tsv"
map   = opt$map # map = "/Users/khaled/Desktop/bioinfo_snakemake/data/map/kfp5s2_CPd_GAPd_HNS_LAPd.txt"
threshold = opt$threshold # threshold = "0.8"
threshold = as.numeric(threshold)
categories = opt$categories # categories = "ShortCond"
output = opt$output # output = "/Users/khaled/Desktop/bioinfo_snakemake/data/network/Func_GAPd_LAPd_CPd_HNS_subsystem/"

# Read the files
myOTU = import_biom(input)

# Standardize the name of the desired category in the mapping file. This makes subsequent operations easier
myMap = import_qiime_sample_data(map)
myMap = data.frame(myMap)
colnames(myMap)[colnames(myMap) == categories] <- "myCat"
myCat = unique(as.vector(myMap$myCat))

for(currentCat in myCat){ # This is what will loop on all categories of your mapping file. Disabled now as the category is being passed on from optparse as the "group" variable
  
  # Filter by condition
  currentSamples = dplyr::filter(myMap, myCat == currentCat) %>% .$X.SampleID %>% as.character(.)
  currentOTU = prune_samples(currentSamples, myOTU) # Filter by condition
  
  # Filter by core threshold
  currentthreshold = length(currentSamples) - (threshold * length(currentSamples))
  
  #There is prune_taxa
  myCurrentOTU = currentOTU@otu_table %>% as.data.frame(.)
  ID = rownames(myCurrentOTU)
  myCurrentOTU = cbind(ID,myCurrentOTU)
  zeroCounts = myCurrentOTU%>%
    gather(var, val, -ID) %>% #Transforming the data from wide to long format
    group_by(val, ID) %>% #Grouping 
    summarise(count = n()) %>% #Performing the count
    reshape2::dcast(ID~val, value.var = "count", fill = 0) #Reshaping the data
  
  colnames(zeroCounts)[colnames(zeroCounts) == "0"] <- "My0"
  
  coreFeatures = zeroCounts %>% select(ID, My0) %>% .[.[, "My0"] < currentthreshold,] %>% .$ID #%>% as.character(.)
  
  YesOrNo = as.character(myCurrentOTU$ID) %in% as.character(coreFeatures)
  
  myTable = cbind(myCurrentOTU,YesOrNo)
  myTable = myTable[myTable$YesOrNo == "TRUE", ]       
  myTable$YesOrNo = NULL # Remove column
  colnames(myTable)[colnames(myTable) == "ID"] <- "#OTU ID" #Change name of column before saving the file to the standard in biom files
  
  
  myOutput = paste0(output,"core+",currentCat,".tsv")
  write_tsv(myTable,myOutput,col_names = TRUE)
}
