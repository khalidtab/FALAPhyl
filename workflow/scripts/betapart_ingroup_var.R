#!/usr/local/bin/Rscript --vanilla



set.seed(1234)
suppressWarnings(suppressMessages(library(optparse)))	
suppressWarnings(suppressMessages(library(betapart)))	
suppressWarnings(suppressMessages(require(tidyverse)))	

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="tsv file", metavar="Features input file formatted as tsv from a biom file"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder to save the files", metavar="Output folder")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

biomFile = opt$input
output = opt$output 

subsetTable = read.delim(biomFile, row.names = 1)

sink(output)
print(paste("# Per Baselga (2012), Bray-Curtis dissimilarity can be broken down to 'Balanced variation in abundance' and 'Abundance gradients'. Simply put, 'balanced variation in abundance' is analogous to 'Spatial turnover' in incidence terms, that is, replacement by another species. 'abundance gradients' means that one sample is a subset of the other. This is analogous to 'nestness' in incidence terms, that is, removal and no replacement by another species. This script will calculate raw Bray-Curtis dissimilarity, as well as the two components of Bray-Curtis"))

brayCalc.core = betapart.core.abund(subsetTable)
brayCalc.multi = beta.multi.abund(brayCalc.core) # multiple site measures (this will give you the overall values)

print(paste0("Within group values: Bray Curtis: ",brayCalc.multi$beta.BRAY," Balanced: ", brayCalc.multi$beta.BRAY.BAL," Gradients: ", brayCalc.multi$beta.BRAY.GRA))
print(paste0("This means that Balanced represents:", (brayCalc.multi$beta.BRAY.BAL/brayCalc.multi$beta.BRAY)*100,"% and Gradients represents ", (brayCalc.multi$beta.BRAY.GRA/brayCalc.multi$beta.BRAY)*100, "%"))
sink()