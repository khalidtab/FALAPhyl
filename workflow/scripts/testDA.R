#!/usr/local/bin/Rscript --vanilla

set.seed(1234)

library("DAtest")
library("tidyverse")
suppressWarnings(suppressMessages(library("optparse")))	

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Biom file", metavar="Features input file formatted as biom"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-c", "--category"), type="character", default=NULL, help="Category", metavar="Name of category column to compare"),
  make_option(c("-s", "--minsample"), type="character", default=NULL, help="Prefiltering number. Minimal number of samples a feature needs to be present in. Otherwise it will be filtered out, and combined as Others", metavar="Min number of samples"),
  make_option(c("-r", "--minread"), type="character", default=NULL, help="Prefiltering number. Minimal number of reads a feature needs to be present in. Otherwise it will be filtered out, and combined as Others", metavar="Min number of reads"),
  make_option(c("-a", "--minabund"), type="character", default=NULL, help="Prefiltering number. Minimal mean relative abundance a feature needs to be present in. Otherwise it will be filtered out, and combined as Others", metavar="Min number of mean relative abundance"),
  make_option(c("-t", "--thetest"), type="character", default=NULL, help="Test to be run", metavar="Test to be run"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file name", metavar="Output file name"),
  make_option(c("-g", "--graph"), type="character", default=NULL, help="Output graph file name", metavar="Output graph file name")
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

mymethod = opt$thetest

if (mymethod == "abc"){
  final=testDA(df, predictor = vec,tests=c("kru","abc"))
  
  
} else if (mymethod == "adx") {
  final=testDA(df, predictor = vec,tests=c("kru","adx"))
  
  
} else if (mymethod == "aov") {
  final=testDA(df, predictor = vec,tests=c("kru","aov"))
  
  
} else if (mymethod == "lao") {
  final=testDA(df, predictor = vec,tests=c("kru","lao"))
  
  
}  else if (mymethod == "lao2") {
  final=testDA(df, predictor = vec,tests=c("kru","lao2"))
  
  
} else if (mymethod == "bay") {
  final=testDA(df, predictor = vec,tests=c("kru","bay"))
  
  
} else if (mymethod == "pea") {
  final=testDA(df, predictor = vec,tests=c("kru","pea"))
  
  
} else if (mymethod == "spe") {
  final=testDA(df, predictor = vec,tests=c("kru","spe"))
  
  
} else if (mymethod == "ds2x") {
  final=testDA(df, predictor = vec,tests=c("kru","ds2x"))
  
  
} else if (mymethod == "ds2") {
  final=testDA(df, predictor = vec,tests=c("kru","ds2"))
  
  
} else if (mymethod == "ere") {
  final=testDA(df, predictor = vec,tests=c("kru","ere"))
  
  
} else if (mymethod == "ere2") {
  final=testDA(df, predictor = vec,tests=c("kru","ere2"))
  
  
} else if (mymethod == "erq") {
  final=testDA(df, predictor = vec,tests=c("kru","erq"))
  
  
} else if (mymethod == "erq2") {
  final=testDA(df, predictor = vec,tests=c("kru","erq2"))
  
  
} else if (mymethod == "fri") {
  final=testDA(df, predictor = vec,tests=c("kru","fri"))
  
  
} else if (mymethod == "neb") {
  final=testDA(df, predictor = vec,tests=c("kru","neb"))
  
  
} else if (mymethod == "poi") {
  final=testDA(df, predictor = vec,tests=c("kru","poi"))
  
  
} else if (mymethod == "qpo") {
  final=testDA(df, predictor = vec,tests=c("kru","qpo"))
  
  
} else if (mymethod == "znb") {
  final=testDA(df, predictor = vec,tests=c("kru","znb"))
  
  
} else if (mymethod == "zpo") {
  final=testDA(df, predictor = vec,tests=c("kru","zpo"))
  
  
} else if (mymethod == "lim") {
  final=testDA(df, predictor = vec,tests=c("kru","lim"))
  
  
} else if (mymethod == "lli") {
  final=testDA(df, predictor = vec,tests=c("kru","lli"))
  
  
} else if (mymethod == "lli2") {
  final=testDA(df, predictor = vec,tests=c("kru","lli2"))
  
  
} else if (mymethod == "vli") {
  final=testDA(df, predictor = vec,tests=c("kru","vli"))
  
  
} else if (mymethod == "lrm") {
  final=testDA(df, predictor = vec,tests=c("kru","lrm"))
  
  
} else if (mymethod == "llm") {
  final=testDA(df, predictor = vec,tests=c("kru","llm"))
  
  
} else if (mymethod == "llm2") {
  final=testDA(df, predictor = vec,tests=c("kru","llm2"))
  
  
} else if (mymethod == "msf") {
  final=testDA(df, predictor = vec,tests=c("kru","msf"))
  
  
} else if (mymethod == "zig") {
  final=testDA(df, predictor = vec,tests=c("kru","zig"))
  
  
} else if (mymethod == "mva") {
  final=testDA(df, predictor = vec,tests=c("kru","mva"))
  
  
} else if (mymethod == "per") {
  final=testDA(df, predictor = vec,tests=c("kru","per"))
  
  
} else if (mymethod == "qua") {
  final=testDA(df, predictor = vec,tests=c("kru","qua"))
  
  
} else if (mymethod == "sam") {
  final=testDA(df, predictor = vec,tests=c("kru","sam"))
  
  
} else if (mymethod == "ttt") {
  final=testDA(df, predictor = vec,tests=c("kru","ttt"))
  
  
} else if (mymethod == "ltt") {
  final=testDA(df, predictor = vec,tests=c("kru","ltt"))
  
  
} else if (mymethod == "ltt2") {
  final=testDA(df, predictor = vec,tests=c("kru","ltt2"))
  
  
} else if (mymethod == "ttr") {
  final=testDA(df, predictor = vec,tests=c("kru","ttr"))
  
  
} else if (mymethod == "wil") {
  final=testDA(df, predictor = vec,tests=c("kru","wil"))
  
  
}  else if (mymethod == "ttc") {
  final=testDA(df, predictor = vec,tests=c("kru","ttc"))
  
  
}  else if (mymethod == "lic") {
  final=testDA(df, predictor = vec,tests=c("kru","lic"))
  
  
} else if (mymethod == "tta") {
  final=testDA(df, predictor = vec,tests=c("kru","tta"))
  
  
} else if (mymethod == "lia") {
  final=testDA(df, predictor = vec,tests=c("kru","lia"))
  
}

write_tsv(as.data.frame(final$table),opt$output)
