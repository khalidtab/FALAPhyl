#!/usr/local/bin/Rscript --vanilla

set.seed(1234)

library("DAtest")
library("tidyverse")
suppressWarnings(suppressMessages(library("optparse")))	

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Biom file", metavar="Features input file formatted as biom"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-t", "--test"), type="character", default=NULL, help="DAtest method", metavar="Type of test that DAtest does"),
  make_option(c("-c", "--category"), type="character", default=NULL, help="Category", metavar="Name of category column to compare"),
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

mymethod = opt$test

if (mymethod == "abc"){
  final=powerDA(df, predictor = vec, test = "abc", cores = 1)
  
  
} else if (mymethod == "adx") {
  final=powerDA(df, predictor = vec, test = "adx", cores = 1)
  
  
} else if (mymethod == "aov") {
  final=powerDA(df, predictor = vec, test = "aov", cores = 1)
  
  
} else if (mymethod == "lao") {
  final=powerDA(df, predictor = vec, test = "lao", cores = 1)
  
  
}  else if (mymethod == "lao2") {
  final=powerDA(df, predictor = vec, test = "lao2", cores = 1)
  
  
} else if (mymethod == "bay") {
  final=powerDA(df, predictor = vec, test = "bay", cores = 1)
  
  
} else if (mymethod == "pea") {
  final=powerDA(df, predictor = vec, test = "pea", cores = 1)
  
  
} else if (mymethod == "spe") {
  final=powerDA(df, predictor = vec, test = "spe", cores = 1)
  
  
} else if (mymethod == "ds2x") {
  final=powerDA(df, predictor = vec, test = "ds2x", cores = 1)
  
  
} else if (mymethod == "ds2") {
  final=powerDA(df, predictor = vec, test = "ds2", cores = 1)
  
  
} else if (mymethod == "ere") {
  final=powerDA(df, predictor = vec, test = "ere", cores = 1)
  
  
} else if (mymethod == "ere2") {
  final=powerDA(df, predictor = vec, test = "ere2", cores = 1)
  
  
} else if (mymethod == "erq") {
  final=powerDA(df, predictor = vec, test = "erq", cores = 1)
  
  
} else if (mymethod == "erq2") {
  final=powerDA(df, predictor = vec, test = "erq2", cores = 1)
  
  
} else if (mymethod == "fri") {
  final=powerDA(df, predictor = vec, test = "fri", cores = 1)
  
  
} else if (mymethod == "neb") {
  final=powerDA(df, predictor = vec, test = "neb", cores = 1)
  
  
} else if (mymethod == "poi") {
  final=powerDA(df, predictor = vec, test = "poi", cores = 1)
  
  
} else if (mymethod == "qpo") {
  final=powerDA(df, predictor = vec, test = "qpo", cores = 1)
  
  
} else if (mymethod == "znb") {
  final=powerDA(df, predictor = vec, test = "znb", cores = 1)
  
  
} else if (mymethod == "zpo") {
  final=powerDA(df, predictor = vec, test = "zpo", cores = 1)
  
  
} else if (mymethod == "kru") {
  final=powerDA(df, predictor = vec, test = "kru", cores = 1)
  
  
} else if (mymethod == "lim") {
  final=powerDA(df, predictor = vec, test = "lim", cores = 1)
  
  
} else if (mymethod == "lli") {
  final=powerDA(df, predictor = vec, test = "lli", cores = 1)
  
  
} else if (mymethod == "lli2") {
  final=powerDA(df, predictor = vec, test = "lli2", cores = 1)
  
  
} else if (mymethod == "vli") {
  final=powerDA(df, predictor = vec, test = "vli", cores = 1)
  
  
} else if (mymethod == "lrm") {
  final=powerDA(df, predictor = vec, test = "lrm", cores = 1)
  
  
} else if (mymethod == "llm") {
  final=powerDA(df, predictor = vec, test = "llm", cores = 1)
  
  
} else if (mymethod == "llm2") {
  final=powerDA(df, predictor = vec, test = "llm2", cores = 1)
  
  
} else if (mymethod == "msf") {
  final=powerDA(df, predictor = vec, test = "msf", cores = 1)
  
  
} else if (mymethod == "zig") {
  final=powerDA(df, predictor = vec, test = "zig", cores = 1)
  
  
} else if (mymethod == "mva") {
  final=powerDA(df, predictor = vec, test = "mva", cores = 1)
  
  
} else if (mymethod == "per") {
  final=powerDA(df, predictor = vec, test = "per", cores = 1)
  
  
} else if (mymethod == "qua") {
  final=powerDA(df, predictor = vec, test = "qua", cores = 1)
  
  
} else if (mymethod == "sam") {
  final=powerDA(df, predictor = vec, test = "sam", cores = 1) 
  
  
} else if (mymethod == "ttt") {
  final=powerDA(df, predictor = vec, test = "ttt", cores = 1)
  
  
} else if (mymethod == "ltt") {
  final=powerDA(df, predictor = vec, test = "ltt", cores = 1)
  
  
} else if (mymethod == "ltt2") {
  final=powerDA(df, predictor = vec, test = "ltt2", cores = 1)
  
  
} else if (mymethod == "ttr") {
  final=powerDA(df, predictor = vec, test = "ttr", cores = 1)
  
  
} else if (mymethod == "wil") {
  final=powerDA(df, predictor = vec, test = "wil", cores = 1)
  
  
}  else if (mymethod == "ttc") {
  final=powerDA(df, predictor = vec, test = "ttc", cores = 1)
  
  
}  else if (mymethod == "lic") {
  final=powerDA(df, predictor = vec, test = "lic", cores = 1)
  
  
 } else if (mymethod == "tta") {
  final=powerDA(df, predictor = vec, test = "tta", cores = 1)
  
  
 } else if (mymethod == "lia") {
   final=powerDA(df, predictor = vec, test = "lia", cores = 1)
   
 }

write_tsv(summary(final),opt$output)
