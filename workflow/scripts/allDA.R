#!/usr/local/bin/Rscript --vanilla

set.seed(1234)

library("DAtest")
library("tidyverse")
suppressWarnings(suppressMessages(library("optparse")))	

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Biom file", metavar="Features input file formatted as biom"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-c", "--category"), type="character", default=NULL, help="Category", metavar="Name of category column to compare"),
  make_option(c("-g", "--graph"), type="character", default=NULL, help="output file name for the graph", metavar="Output file name for the graph"),
  make_option(c("-s", "--minsample"), type="character", default=NULL, help="Prefiltering number. Minimal number of samples a feature needs to be present in. Otherwise it will be filtered out, and combined as Others", metavar="Min number of samples"),
  make_option(c("-r", "--minread"), type="character", default=NULL, help="Prefiltering number. Minimal number of reads a feature needs to be present in. Otherwise it will be filtered out, and combined as Others", metavar="Min number of reads"),
  make_option(c("-a", "--minabund"), type="character", default=NULL, help="Prefiltering number. Minimal mean relative abundance a feature needs to be present in. Otherwise it will be filtered out, and combined as Others", metavar="Min number of mean relative abundance"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file name", metavar="Output file name")
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
if ( as.numeric(opt$minsample)==0 && as.numeric(opt$minread)==0 && as.numeric(opt$minabund)==0){df = df}else{
  df = preDA(df, min.samples = as.numeric(opt$minsample), min.reads = as.numeric(opt$minread), min.abundance =  as.numeric(opt$minabund))
}
# Load mapping file
map = opt$mapping
map = read.csv(map,sep="\t") %>% as.data.frame(.)

category = opt$category
catNum = which(colnames(map) == category)
working_map = cbind(as.character(map[,1]),
                    as.character(map[,catNum])) %>% as.data.frame(.)
colnames(working_map) = c("SampleID","condition")
vec = working_map$condition %>% as.factor(.)

testa = allDA(df, predictor = vec, cores = 2)

output = opt$output

for (x in 1:length(testa$results)){
write_delim(as.data.frame(testa$results[x]),file=output,delim="\t",append=TRUE,col_names=TRUE)
}
graphOutput = opt$graph

svg(filename = graphOutput,family="sans")
plota
dev.off()
