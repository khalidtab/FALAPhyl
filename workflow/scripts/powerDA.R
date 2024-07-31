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
  make_option(c("-l", "--log"), type="character", default=NULL, help="log file", metavar="log file"), 
   make_option(c("-e", "--tmp"), type="character", default=NULL, help="tmp folder", metavar="tmp folder"), 
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

tmp_folder = opt$tmp

log_file = opt$log

mymethod = opt$test

if (mymethod == "abc"){
  final=powerDA(df, predictor = vec, test = "abc")
  
  
} else if (mymethod == "adx") {
  final=powerDA(df, predictor = vec, test = "adx")
  
  
} else if (mymethod == "aov") {
  final=powerDA(df, predictor = vec, test = "aov")
  
  
} else if (mymethod == "lao") {
  final=powerDA(df, predictor = vec, test = "lao")
  
  
}  else if (mymethod == "lao2") {
  final=powerDA(df, predictor = vec, test = "lao2")
  
  
} else if (mymethod == "bay") {
  final=powerDA(df, predictor = vec, test = "bay")
  
  
} else if (mymethod == "pea") {
  final=powerDA(df, predictor = vec, test = "pea")
  
  
} else if (mymethod == "spe") {
  final=powerDA(df, predictor = vec, test = "spe")
  
  
} else if (mymethod == "ds2x") {
  final=powerDA(df, predictor = vec, test = "ds2x")
  
  
} else if (mymethod == "ds2") {
  final=powerDA(df, predictor = vec, test = "ds2")
  
  
} else if (mymethod == "ere") {
  final=powerDA(df, predictor = vec, test = "ere")
  
  
} else if (mymethod == "ere2") {
  final=powerDA(df, predictor = vec, test = "ere2")
  
  
} else if (mymethod == "erq") {
  final=powerDA(df, predictor = vec, test = "erq")
  
  
} else if (mymethod == "erq2") {
  final=powerDA(df, predictor = vec, test = "erq2")
  
  
} else if (mymethod == "fri") {
  final=powerDA(df, predictor = vec, test = "fri")
  
  
} else if (mymethod == "neb") {
  final=powerDA(df, predictor = vec, test = "neb")
  
  
} else if (mymethod == "poi") {
  final=powerDA(df, predictor = vec, test = "poi")
  
  
} else if (mymethod == "qpo") {
  final=powerDA(df, predictor = vec, test = "qpo")
  
  
} else if (mymethod == "znb") {
  final=powerDA(df, predictor = vec, test = "znb")
  
  
} else if (mymethod == "zpo") {
  final=powerDA(df, predictor = vec, test = "zpo")
  
  
} else if (mymethod == "kru") {
  final=powerDA(df, predictor = vec, test = "kru")
  
  
} else if (mymethod == "lim") {
  final=powerDA(df, predictor = vec, test = "lim")
  
  
} else if (mymethod == "lli") {
  final=powerDA(df, predictor = vec, test = "lli")
  
  
} else if (mymethod == "lli2") {
  final=powerDA(df, predictor = vec, test = "lli2")
  
  
} else if (mymethod == "vli") {
  final=powerDA(df, predictor = vec, test = "vli")
  
  
} else if (mymethod == "lrm") {
  final=powerDA(df, predictor = vec, test = "lrm")
  
  
} else if (mymethod == "llm") {
  final=powerDA(df, predictor = vec, test = "llm")
  
  
} else if (mymethod == "llm2") {
  final=powerDA(df, predictor = vec, test = "llm2")
  
  
} else if (mymethod == "msf") {
  final=powerDA(df, predictor = vec, test = "msf")
  
  
} else if (mymethod == "zig") {
  final=powerDA(df, predictor = vec, test = "zig")
  
  
} else if (mymethod == "mva") {
  final=powerDA(df, predictor = vec, test = "mva")
  
  
} else if (mymethod == "per") {
  final=powerDA(df, predictor = vec, test = "per")
  
  
} else if (mymethod == "qua") {
  final=powerDA(df, predictor = vec, test = "qua")
  
  
} else if (mymethod == "sam") {
  final=powerDA(df, predictor = vec, test = "sam") 
  
  
} else if (mymethod == "ttt") {
  final=powerDA(df, predictor = vec, test = "ttt")
  
  
} else if (mymethod == "ltt") {
  final=powerDA(df, predictor = vec, test = "ltt")
  
  
} else if (mymethod == "ltt2") {
  final=powerDA(df, predictor = vec, test = "ltt2")
  
  
} else if (mymethod == "ttr") {
  final=powerDA(df, predictor = vec, test = "ttr")
  
  
} else if (mymethod == "wil") {
  final=powerDA(df, predictor = vec, test = "wil")
  
  
}  else if (mymethod == "ttc") {
  final=powerDA(df, predictor = vec, test = "ttc")
  
  
}  else if (mymethod == "lic") {
  final=powerDA(df, predictor = vec, test = "lic")
  
  
 } else if (mymethod == "tta") {
  final=powerDA(df, predictor = vec, test = "tta")
  
  
 } else if (mymethod == "lia") {
   final=powerDA(df, predictor = vec, test = "lia")
   
 } else if (mymethod %in% c("CPLM","ZICP","ZSCP","ZACP")){
   
   final = NA
   final$text = "effect size power calculation is not available for the chosen test"
   
 }

write_tsv(summary(final),opt$output)
