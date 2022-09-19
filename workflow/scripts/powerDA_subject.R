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
  make_option(c("-p", "--pair"), type="character", default=NULL, help="Column name of the subject ID in the mapping file", metavar="Column name of the subject ID in the mapping file"),  
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
subjectID = opt$pair

catNum = which(colnames(map) == category)
subjectNum = which(colnames(map) == subjectID)
working_map = cbind(as.character(map[,1]),
                    as.character(map[,catNum]),
                    as.character(map[,subjectNum])) %>% as.data.frame(.)
colnames(working_map) = c("SampleID","condition","subject")
vec = working_map$condition %>% as.factor(.)
subject = working_map$subject %>% as.factor(.)

mymethod = opt$test

if (mymethod == "fri"){final=powerDA(df, predictor = vec, test = "fri", paired=subject)
} else if (mymethod == "qua"){final=powerDA(df, predictor = vec, test = "qua", paired=subject)
} else if (mymethod == "per"){final=powerDA(df, predictor = vec, test = "per", paired=subject)
} else if (mymethod == "sam"){final=powerDA(df, predictor = vec, test = "sam", paired=subject)
} else if (mymethod == "ttt"){final=powerDA(df, predictor = vec, test = "ttt", paired=subject)
} else if (mymethod == "ttt"){final=powerDA(df, predictor = vec, test = "ttt", paired=subject)
} else if (mymethod == "ttt"){final=powerDA(df, predictor = vec, test = "ttt", paired=subject)
} else if (mymethod == "ltt"){final=powerDA(df, predictor = vec, test = "ltt", paired=subject)
} else if (mymethod == "ltt2"){final=powerDA(df, predictor = vec, test = "ltt2", paired=subject)
} else if (mymethod == "ttr"){final=powerDA(df, predictor = vec, test = "ttr", paired=subject)
} else if (mymethod == "wil"){final=powerDA(df, predictor = vec, test = "wil", paired=subject)
} else if (mymethod == "ds2x"){final=powerDA(df, predictor = vec, test = "ds2x", paired=subject)
} else if (mymethod == "ds2"){final=powerDA(df, predictor = vec, test = "ds2", paired=subject)
} else if (mymethod == "erq"){final=powerDA(df, predictor = vec, test = "erq", paired=subject)
} else if (mymethod == "erq2"){final=powerDA(df, predictor = vec, test = "erq2", paired=subject)
} else if (mymethod == "mva"){final=powerDA(df, predictor = vec, test = "mva", paired=subject)
} else if (mymethod == "neb"){final=powerDA(df, predictor = vec, test = "neb", paired=subject)
} else if (mymethod == "poi"){final=powerDA(df, predictor = vec, test = "poi", paired=subject)
} else if (mymethod == "lim"){final=powerDA(df, predictor = vec, test = "lim", paired=subject)
} else if (mymethod == "lim"){final=powerDA(df, predictor = vec, test = "lim", paired=subject)
} else if (mymethod == "lim"){final=powerDA(df, predictor = vec, test = "lim", paired=subject)
} else if (mymethod == "lli"){final=powerDA(df, predictor = vec, test = "lli", paired=subject)
} else if (mymethod == "lli2"){final=powerDA(df, predictor = vec, test = "lli2", paired=subject)
} else if (mymethod == "vli"){final=powerDA(df, predictor = vec, test = "vli", paired=subject)
} else if (mymethod == "lrm"){final=powerDA(df, predictor = vec, test = "lrm", paired=subject)
} else if (mymethod == "lrm"){final=powerDA(df, predictor = vec, test = "lrm", paired=subject)
} else if (mymethod == "lrm"){final=powerDA(df, predictor = vec, test = "lrm", paired=subject)
} else if (mymethod == "llm"){final=powerDA(df, predictor = vec, test = "llm", paired=subject)
} else if (mymethod == "llm2"){final=powerDA(df, predictor = vec, test = "llm2", paired=subject)
} else if (mymethod == "zig"){final=powerDA(df, predictor = vec, test = "zig", paired=subject)}

write_tsv(summary(final),opt$output)
