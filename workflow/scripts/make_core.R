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


for(currentCat in myCat){
currentSamples = dplyr::filter(myMap, myCat == currentCat) %>% .$X.SampleID %>% as.character(.)
myCurrentOTU = myOTU[ , names(myOTU) %in% currentSamples] 
myCurrentOTU[myCurrentOTU==0] = NA #Makes it easier to work on the samples
myCurrentOTU = myCurrentOTU %>% filter(if_any(everything(), ~ !is.na(.))) # Rows with only NAs are removed


retainedOTUs = myCurrentOTU[which(rowMeans(!is.na(myCurrentOTU)) >= threshold), ] # Convert matrix to TRUE and FALSE. True is represented by 1 and false by 0. Get the mean of the row. If that row is more than the threshold then it is retained

if(dim(retainedOTUs)[1] == 0){ # ie, no retained OTUs because of the threshold
  print(paste0("The level of core you have specified resulted in no retained features for",currentCat,". Adjust your level of core and rerun snakemake on a lower core level 'threshold' by adjusting the input.yaml file."))
}

excludedOTUs = subset(myCurrentOTU, !row.names(myCurrentOTU) %in% row.names(retainedOTUs)) # Give me myOTUs were row names of myOTUs is NOT (!) in row names of retainedOTUS
excludedOTUs = colSums(excludedOTUs,na.rm = TRUE) %>% as.matrix(.) %>% t(.) # The sum of all OTUs that did not pass the filtering step. This is done to reduce zeros in the matrix, and to maintain the data as it is compositional in nature. After SparCC, it can be deleted due to to satisfy the subcompositional coherence
row.names(excludedOTUs) = "Others"
excludedOTUs = as.data.frame(excludedOTUs)
myCurrentOTU_amalgamated = rbind(retainedOTUs,excludedOTUs)



ID = rownames(myCurrentOTU_amalgamated)
myCurrentOTU_amalgamated = cbind(ID,myCurrentOTU_amalgamated)
colnames(myCurrentOTU_amalgamated)[colnames(myCurrentOTU_amalgamated) == "ID"] <- "#OTU ID" #Change name of column before saving the file to the standard in biom files
myCurrentOTU_amalgamated[is.na(myCurrentOTU_amalgamated)] = 0

myOutput = paste0(output,currentCat,".tsv")
write_tsv(myCurrentOTU_amalgamated,myOutput,col_names = TRUE)
}
