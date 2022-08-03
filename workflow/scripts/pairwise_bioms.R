#!/usr/local/bin/Rscript --vanilla



set.seed(1234)
suppressWarnings(suppressMessages(library(optparse)))	
suppressWarnings(suppressMessages(require(tidyverse)))	
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(phyloseq)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Biom file", metavar="Features input file formatted as biom"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-g", "--groups"), type="character", default=NULL, help="Column name in the mapping file that indicates which grouping the samples belong to", metavar="Group indication column"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder to save the files", metavar="Output folder")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

biomFile = opt$input 
mappingFile = opt$mapping
condition = opt$groups 
output = opt$output 

####### Functions #######
load_phylo = function(myPhyloseqObject){
  # Load and format tables
  myBiomTSV = myPhyloseqObject@otu_table@.Data
  
  featureNames = row.names(myBiomTSV) %>% as.matrix()
  myBiomTSV = as.data.frame(myBiomTSV)
  row.names(myBiomTSV) = featureNames
  myBiomTSV = myBiomTSV[-c(1)] %>% t(.) %>% as.data.frame(.)
  return(myBiomTSV)
}
#########################

# Make the phyloseq merged file


biom = suppressMessages(readr::read_tsv(biomFile,skip=1)) %>% as.data.frame(.)
rownames(biom) = biom[,1]
biom[,1] = NULL

biom = otu_table(biom, taxa_are_rows = TRUE)

# Now do the required calculations on each group at a time
## First, find what the groups are
mapping = import_qiime_sample_data(mappingFile)
condNum = which(colnames(mapping)==condition)
condName = mapping[,condNum][[1]] %>% unique(.) %>% as.character(.) %>% combn(.,2,simplify=FALSE)

myPhylo = merge_phyloseq(biom,mapping)
myCondPop = c()

for (myCond in condName) {
subsetPhylo = subset_samples(myPhylo, myPhylo@sam_data@.Data[[condNum]]==myCond[1])
subsetPhylo2 = subset_samples(myPhylo, myPhylo@sam_data@.Data[[condNum]]==myCond[2])

subsetTable = suppressWarnings(suppressMessages(load_phylo(subsetPhylo)))
subsetTable2 = suppressWarnings(suppressMessages(load_phylo(subsetPhylo2)))

subsetMeta  = sample_data(subsetPhylo)  %>% as.matrix(.) %>% as.data.frame(.)
subsetMeta2 = sample_data(subsetPhylo2) %>% as.matrix(.) %>% as.data.frame(.)

combinedBiom = bind_rows(subsetTable, subsetTable2)
combinedmeta = bind_rows(subsetMeta,subsetMeta2)

write.table(combinedBiom,file=paste0(output,"/",condition,"–",myCond[1],"+",myCond[2],".tsv"), 
            quote = FALSE, sep = "\t",col.names = NA)

write.table(combinedmeta,file=paste0(output,"/",condition,"–",myCond[1],"+",myCond[2],".txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
}
