#!/usr/local/bin/Rscript --vanilla

suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(dunn.test)))
suppressWarnings(suppressMessages(library(magrittr)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="AUC/FDR/Power text", metavar="AUC/FDR/Power text"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder for all pairwise files", metavar="Output folder name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

myFile = read_tsv(opt$input)

if (length(unique(myFile$Method)) == 2) {
  myResults = wilcox.test(score ~ Method, data = myFile)
} else if (length(unique(myFile$Method)) > 2) {
  myResults = dunn.test(myFile$score, myFile$Method, method = "bh")
} else {
  message("Only one group found. No test needed.")
  myResults = as.data.frame("one group")
}

write_tsv(myResults,opt$output)