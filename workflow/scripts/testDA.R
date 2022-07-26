#!/usr/local/bin/Rscript --vanilla

set.seed(1234)

library("DAtest")
library("tidyverse")
suppressWarnings(suppressMessages(library("optparse")))	

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Biom file", metavar="Features input file formatted as biom"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-c", "--category"), type="character", default=NULL, help="Category", metavar="Name of category column to compare"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file name", metavar="Output file name"),
  make_option(c("-g", "--graph"), type="character", default=NULL, help="output file name for the graph", metavar="Output graph file name in svg")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load the biom tsv file
df = opt$input
df = read_tsv(df)
dfRows = as.data.frame(df[,1])
df = as.data.frame(df)
rownames(df) = dfRows[,1]
df[,1] = NULL
df[] = lapply(df, as.numeric)
df = t(df)
# Load mapping file
map = opt$mapping
map = read.csv(map,sep="\t") %>% as.data.frame(.)

category = opt$category
catNum = which(colnames(map) == category)
working_map = cbind(as.character(map[,1]),
                    as.character(map[,catNum])) %>% as.data.frame(.)
colnames(working_map) = c("SampleID","condition")
vec = working_map$condition %>% as.factor(.)

test = testDA(df, predictor = vec)
summary(test)

write.csv(summary(test),output)
svg(opt$graph)
plot(test)
dev.off()
