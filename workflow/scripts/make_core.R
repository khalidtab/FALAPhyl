suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(phyloseq)))
suppressWarnings(suppressMessages(library(tidyverse)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input biom file", metavar="Input biom file"),
  make_option(c("-m", "--map"), type="character", default=NULL, help="The mapping file", metavar="Mapping file"),
  make_option(c("-c", "--categories"), type="character", default=NULL, help="The categories to filter by", metavar="Categories"),
  make_option(c("-t", "--threshold"), type="character", default=NULL, help="Threshold for core filtering", metavar="Threshold"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder for core files", metavar="output")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

input = opt$input 
map   = opt$map #
threshold = opt$threshold 
threshold = as.numeric(threshold)
categories = opt$categories
output = opt$output 

# Read the files
myOTU = import_biom(input) %>% as.data.frame(.)

# Standardize the name of the desired category in the mapping file. This makes subsequent operations easier
myMap = import_qiime_sample_data(map)
myMap = data.frame(myMap)
colnames(myMap)[colnames(myMap) == categories] <- "myCat"
myCat = unique(as.vector(myMap$myCat))

for(currentCat in myCat){ # This is what will loop on all categories of your mapping file. Disabled now as the category is being passed on from optparse as the "group" variable
  
  # Filter by condition
  currentSamples = dplyr::filter(myMap, myCat == currentCat) %>% .$X.SampleID %>% as.character(.)
  myCurrentOTU = myOTU[ , names(myOTU) %in% currentSamples]
  
  #There is prune_taxa
  myCurrentOTU[myCurrentOTU==0] = NA #Makes it easier to work on the samples
  myCurrentOTU = myCurrentOTU %>% filter(if_any(everything(), ~ !is.na(.))) # Rows with only NAs are removed
  myTable = myCurrentOTU[which(rowMeans(!is.na(myCurrentOTU)) > threshold), ]
  
  if(dim(myTable)[1] == 0){ # ie, no retained OTUs because of the threshold
    print(paste0("The level of core you have specified resulted in no retained features for",currentCat,". Adjust your level of core and rerun snakemake on a lower core level 'threshold' by adjusting the input.yaml file."))
  }
  
  ID = rownames(myTable)
  myTable = cbind(ID,myTable)
  colnames(myTable)[colnames(myTable) == "ID"] <- "#OTU ID" #Change name of column before saving the file to the standard in biom files
  myTable[is.na(myTable)] = 0
  
  myOutput = paste0(output,"_core+",currentCat,".tsv")
  write_tsv(myTable,myOutput,col_names = TRUE)
}
