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
  final=testDA(df, predictor = vec,tests=c("kru","abc"),cores=1)
  
  
} else if (mymethod == "adx") {
  final=testDA(df, predictor = vec,tests=c("kru","adx"),cores=1)
  
  
} else if (mymethod == "aov") {
  final=testDA(df, predictor = vec,tests=c("kru","aov"),cores=1)
  
  
} else if (mymethod == "lao") {
  final=testDA(df, predictor = vec,tests=c("kru","lao"),cores=1)
  
  
}  else if (mymethod == "lao2") {
  final=testDA(df, predictor = vec,tests=c("kru","lao2"),cores=1)
  
  
} else if (mymethod == "bay") {
  final=testDA(df, predictor = vec,tests=c("kru","bay"),cores=1)
  
  
} else if (mymethod == "pea") {
  final=testDA(df, predictor = vec,tests=c("kru","pea"),cores=1)
  
  
} else if (mymethod == "spe") {
  final=testDA(df, predictor = vec,tests=c("kru","spe"),cores=1)
  
  
} else if (mymethod == "ds2x") {
  final=testDA(df, predictor = vec,tests=c("kru","ds2x"),cores=1)
  
  
} else if (mymethod == "ds2") {
  final=testDA(df, predictor = vec,tests=c("kru","ds2"),cores=1)
  
  
} else if (mymethod == "ere") {
  final=testDA(df, predictor = vec,tests=c("kru","ere"),cores=1)
  
  
} else if (mymethod == "ere2") {
  final=testDA(df, predictor = vec,tests=c("kru","ere2"),cores=1)
  
  
} else if (mymethod == "erq") {
  final=testDA(df, predictor = vec,tests=c("kru","erq"),cores=1)
  
  
} else if (mymethod == "erq2") {
  final=testDA(df, predictor = vec,tests=c("kru","erq2"),cores=1)
  
  
} else if (mymethod == "fri") {
  final=testDA(df, predictor = vec,tests=c("kru","fri"),cores=1)
  
  
} else if (mymethod == "neb") {
  final=testDA(df, predictor = vec,tests=c("kru","neb"),cores=1)
  
  
} else if (mymethod == "poi") {
  final=testDA(df, predictor = vec,tests=c("kru","poi"),cores=1)
  
  
} else if (mymethod == "qpo") {
  final=testDA(df, predictor = vec,tests=c("kru","qpo"),cores=1)
  
  
} else if (mymethod == "znb") {
  final=testDA(df, predictor = vec,tests=c("kru","znb"),cores=1)
  
  
} else if (mymethod == "zpo") {
  final=testDA(df, predictor = vec,tests=c("kru","zpo"),cores=1)
  
  
} else if (mymethod == "lim") {
  final=testDA(df, predictor = vec,tests=c("kru","lim"),cores=1)
  
  
} else if (mymethod == "lli") {
  final=testDA(df, predictor = vec,tests=c("kru","lli"),cores=1)
  
  
} else if (mymethod == "lli2") {
  final=testDA(df, predictor = vec,tests=c("kru","lli2"),cores=1)
  
  
} else if (mymethod == "vli") {
  final=testDA(df, predictor = vec,tests=c("kru","vli"),cores=1)
  
  
} else if (mymethod == "lrm") {
  final=testDA(df, predictor = vec,tests=c("kru","lrm"),cores=1)
  
  
} else if (mymethod == "llm") {
  final=testDA(df, predictor = vec,tests=c("kru","llm"),cores=1)
  
  
} else if (mymethod == "llm2") {
  final=testDA(df, predictor = vec,tests=c("kru","llm2"),cores=1)
  
  
} else if (mymethod == "msf") {
  final=testDA(df, predictor = vec,tests=c("kru","msf"),cores=1)
  
  
} else if (mymethod == "zig") {
  final=testDA(df, predictor = vec,tests=c("kru","zig"),cores=1)
  
  
} else if (mymethod == "mva") {
  final=testDA(df, predictor = vec,tests=c("kru","mva"),cores=1)
  
  
} else if (mymethod == "per") {
  final=testDA(df, predictor = vec,tests=c("kru","per"),cores=1)
  
  
} else if (mymethod == "qua") {
  final=testDA(df, predictor = vec,tests=c("kru","qua"),cores=1)
  
  
} else if (mymethod == "sam") {
  final=testDA(df, predictor = vec,tests=c("kru","sam"),cores=1)
  
  
} else if (mymethod == "ttt") {
  final=testDA(df, predictor = vec,tests=c("kru","ttt"),cores=1)
  
  
} else if (mymethod == "ltt") {
  final=testDA(df, predictor = vec,tests=c("kru","ltt"),cores=1)
  
  
} else if (mymethod == "ltt2") {
  final=testDA(df, predictor = vec,tests=c("kru","ltt2"),cores=1)
  
  
} else if (mymethod == "ttr") {
  final=testDA(df, predictor = vec,tests=c("kru","ttr"),cores=1)
  
  
} else if (mymethod == "wil") {
  final=testDA(df, predictor = vec,tests=c("kru","wil"),cores=1)
  
  
}  else if (mymethod == "ttc") {
  final=testDA(df, predictor = vec,tests=c("kru","ttc"),cores=1)
  
  
}  else if (mymethod == "lic") {
  final=testDA(df, predictor = vec,tests=c("kru","lic"),cores=1)
  
  
} else if (mymethod == "tta") {
  final=testDA(df, predictor = vec,tests=c("kru","tta"),cores=1)
  
  
} else if (mymethod == "lia") {
  final=testDA(df, predictor = vec,tests=c("kru","lia"),cores=1)
  
}

write_tsv(as.data.frame(final$table),opt$output)

