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



myOTU[myOTU==0] = NA #Makes it easier to work on the samples
myOTU = myOTU %>% filter(if_any(everything(), ~ !is.na(.))) # Rows with only NAs are removed
retainedOTUs = myOTU[which(rowMeans(!is.na(myOTU)) > threshold), ] # Convert matrix to TRUE and FALSE. True is represented by 1 and false by 0. Get the mean of the row. If that row is more than the threshold then it is retained


if(dim(retainedOTUs)[1] == 0){ # ie, no retained OTUs because of the threshold
  print(paste0("The level of core you have specified resulted in no retained features for",currentCat,". Adjust your level of core and rerun snakemake on a lower core level 'threshold' by adjusting the input.yaml file."))
}
excludedOTUs = subset(myOTU, !row.names(myOTU) %in% row.names(retainedOTUs))
Others = colSums(excludedOTUs,na.rm = TRUE) %>% as.matrix(.) %>% t(.) # The sum of all OTUs that did not pass the filtering step. This is done to reduce zeros in the matrix, and to maintain the data as it is compositional in nature. After SparCC, it can be deleted due to to satisfy the subcompositional coherence
row.names(Others) = "Others"

retainedAndOthers = rbind(retainedOTUs,Others)

for(currentCat in myCat)
  currentSamples = dplyr::filter(myMap, myCat == currentCat) %>% .$X.SampleID %>% as.character(.)
myCurrentOTU = retainedAndOthers[ , names(retainedAndOthers) %in% currentSamples] 

ID = rownames(myCurrentOTU)
myCurrentOTU = cbind(ID,myCurrentOTU)
colnames(myCurrentOTU)[colnames(myCurrentOTU) == "ID"] <- "#OTU ID" #Change name of column before saving the file to the standard in biom files
retainedAndOthers[is.na(retainedAndOthers)] = 0

myOutput = paste0(output,"_core+",currentCat,".tsv")
write_tsv(myTable,myOutput,col_names = TRUE)
}
