#install.packages("workflow/scripts/vegan_2.5-6.tar", repos = NULL, type="source", INSTALL_opts = '--no-lock')
library("optparse")
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="path to tsv file of the biom file", metavar="Input biom in tsv format"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="path to output filtered file", metavar="path to output filtered file"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="The mapping file to filter against", metavar="Mapping file")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

tsv = opt$input
mapping_file = opt$mapping 
output = opt$output


myTSV = read_tsv(tsv, skip = 1) %>% as.data.frame(.)
row.names(myTSV) = myTSV[,1]
myTSV[,1] = NULL

mapping_file = read_tsv(mapping_file)
sampleIDs = mapping_file[,1] %>% unlist(.) %>% as.character(.)

myTSV = myTSV %>% .[,which(colnames(.) %in% sampleIDs)] %>% filter(rowSums(.) != 0)  # Keep samples that are in the mapping file, then keep features which the sum of the rows is not zero
myTSV = cbind(rownames(myTSV),myTSV) # Put the rows names back as a column
colnames(myTSV)[1] = "#OTU ID" 



write_tsv(myTSV,output)

