#!/usr/local/bin/Rscript --vanilla

set.seed(1234)
suppressWarnings(suppressMessages(library(optparse)))	
suppressWarnings(suppressMessages(library(betapart)))	
suppressWarnings(suppressMessages(require(tidyverse)))	
suppressWarnings(suppressMessages(require(vegan)))	
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(phyloseq)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Biom file", metavar="Features input file formatted as biom"),
  make_option(c("-b", "--bray"), type="character", default=NULL, help="output file to save the Bray-Curtis tsv matrix", metavar="Output Bray-Curtis"),
  make_option(c("-r", "--replace"), type="character", default=NULL, help="output file to save the turnover/replacement tsv matrix", metavar="Output replacement"),
  make_option(c("-n", "--noreplace"), type="character", default=NULL, help="output file to save the netedness/no replacement tsv matrix", metavar="Output no replacement")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

biomFile = opt$input
replace = opt$replace
noreplace = opt$noreplace
bray = opt$bray

####### Functions #######
brayBreakdown = function(biomMatrix){
  biomMatrix = biomMatrix %>% t() %>% as.matrix()
  BrayBreakdown = beta.pair.abund(biomMatrix, index.family="bray")
  
  BrayReplacement = BrayBreakdown$beta.bray.bal %>% as.matrix() %>% as.data.frame() 
  BrayNoReplacement = BrayBreakdown$beta.bray.gra %>% as.matrix %>% as.data.frame() 
  BrayCurtis = BrayBreakdown$beta.bray %>% as.matrix %>% as.data.frame() 
  
  myBrays = list(BrayCurtis,BrayReplacement,BrayNoReplacement)
  
  return(myBrays)
  
}
#########################

# Make the phyloseq merged file
biom = suppressWarnings(suppressMessages(import_biom(biomFile)))

# Write the Bray Curtis dissimilarity, balanced and gradients matrices
myBrays = brayBreakdown(biom@otu_table@.Data) 

write.table(x=myBrays[[1]],file=paste0(bray), quote = FALSE, sep = "\t",col.names = NA)
write.table(x=myBrays[[2]],file=paste0(replace), quote = FALSE, sep = "\t",col.names = NA)
write.table(x=myBrays[[3]],file=paste0(noreplace), quote = FALSE, sep = "\t",col.names = NA)