suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(require(tidyverse)))
suppressWarnings(suppressMessages(library(phyloseq)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="biom file in tsv format", metavar=" biom file in tsv format"),
  make_option(c("-a", "--alpha"), type="character", default=NULL, help="Name of the desired alpha diversity. Options: Observed, Chao, ACE, Shannon, Simpson, InvSimpson, Fisher", metavar="Output File name"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output file", metavar="Output File name")
);

biom = opt$input
alpha = opt$alpha
output = opt$output


biomFile = read.table(biom,comment.char = "@",skip=1,sep = "\t",header=TRUE)
row.names(biomFile) = biomFile[,1]
biomFile[,1] = NULL 
biomFile = otu_table(biomFile,taxa_are_rows = TRUE)

myAlpha = phyloseq::estimate_richness(biomFile,split=TRUE,measures=alpha)
myAlpha = cbind(row.names(myAlpha),myAlpha)
colnames(myAlpha)[1] = "Samples"

write_delim(myAlpha,output,delim="\t")