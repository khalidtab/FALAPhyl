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
  make_option(c("-p", "--pair"), type="character", default=NULL, help="Column name of the subject ID in the mapping file", metavar="Column name of the subject ID in the mapping file"),  
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
subjectID = opt$pair

catNum = which(colnames(map) == category)
subjectNum = which(colnames(map) == subjectID)
working_map = cbind(as.character(map[,1]),
                    as.character(map[,catNum]),
                    as.character(map[,subjectNum])) %>% as.data.frame(.)
colnames(working_map) = c("SampleID","condition","subject")
vec = working_map$condition %>% as.factor(.)
subject = working_map$subject %>% as.factor(.)

mymethod = opt$thetest


if (mymethod == "fri"){final=testDA(df, predictor = vec, tests=c("wil","fri"), paired=subject)
} else if (mymethod == "qua"){final=testDA(df, predictor = vec, tests=c("wil","qua"), paired=subject)
} else if (mymethod == "per"){final=testDA(df, predictor = vec, tests=c("wil","per"), paired=subject)
} else if (mymethod == "sam"){final=testDA(df, predictor = vec, tests=c("wil","sam"), paired=subject)
} else if (mymethod == "ttt"){final=testDA(df, predictor = vec, tests=c("wil","ttt"), paired=subject)
} else if (mymethod == "ltt"){final=testDA(df, predictor = vec, tests=c("wil","ltt"), paired=subject)
} else if (mymethod == "ltt2"){final=testDA(df, predictor = vec, tests=c("wil","ltt2"), paired=subject)
} else if (mymethod == "ttr"){final=testDA(df, predictor = vec, tests=c("wil","ttr"), paired=subject)
} else if (mymethod == "wil"){final=testDA(df, predictor = vec, tests=c("wil","wil"), paired=subject)
} else if (mymethod == "ds2x"){final=testDA(df, predictor = vec, tests=c("wil","ds2x"), paired=subject)
} else if (mymethod == "ds2"){final=testDA(df, predictor = vec, tests=c("wil","ds2"), paired=subject)
} else if (mymethod == "erq"){final=testDA(df, predictor = vec, tests=c("wil","erq"), paired=subject)
} else if (mymethod == "erq2"){final=testDA(df, predictor = vec, tests=c("wil","erq2"), paired=subject)
} else if (mymethod == "neb"){final=testDA(df, predictor = vec, tests=c("wil","neb"), paired=subject)
} else if (mymethod == "poi"){final=testDA(df, predictor = vec, tests=c("wil","poi"), paired=subject)
} else if (mymethod == "lim"){final=testDA(df, predictor = vec, tests=c("wil","lim"), paired=subject)
} else if (mymethod == "lli"){final=testDA(df, predictor = vec, tests=c("wil","lli"), paired=subject)
} else if (mymethod == "lli2"){final=testDA(df, predictor = vec, tests=c("wil","lli2"), paired=subject)
} else if (mymethod == "vli"){final=testDA(df, predictor = vec, tests=c("wil","vli"), paired=subject)
} else if (mymethod == "lrm"){final=testDA(df, predictor = vec, tests=c("wil","lrm"), paired=subject)
} else if (mymethod == "llm"){final=testDA(df, predictor = vec, tests=c("wil","llm"), paired=subject)
} else if (mymethod == "llm2"){final=testDA(df, predictor = vec, tests=c("wil","llm2"), paired=subject)
} else if (mymethod == "zig"){final=testDA(df, predictor = vec, tests=c("wil","zig"), paired=subject)
} else if (mymethod %in% c("CPLM","ZICP","ZSCP","ZACP")){
  
  library("Tweedieverse")
  
  df <- apply(df, 2, function(x) x / sum(x))
  
  CPLM <- function(count_table, predictor, paired, covars) {
    
    write_log <- function(message) {
      write(message, file = log_file, append = TRUE)
    }
    
    # Create a unique output directory based on timestamp
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- paste0(tmp_folder, timestamp)
    dir.create(output_dir, recursive = TRUE)
    results_file <- paste0(output_dir, "/all_results.tsv")
    
    write_log("Entering function")
    
    # Debug: Print inputs
    write_log("count_table:")
    write_log(paste(capture.output(head(count_table)), collapse = "\n"))
    write_log(paste("Dimensions of count_table:", paste(dim(count_table), collapse = " x ")))
    
    write_log("predictor:")
    write_log(paste(capture.output(head(predictor)), collapse = "\n"))
    write_log(paste("Length of predictor:", length(predictor)))
    
    # Prepare the predictor data frame
    predictor_df <- data.frame(sampleid = colnames(count_table), metadata = as.character(predictor))
    rownames(predictor_df) <- predictor_df$sampleid
    predictor_df$sampleid <- NULL
    
    write_log("predictor_df:")
    write_log(paste(capture.output(head(predictor_df)), collapse = "\n"))
    write_log(paste("Dimensions of predictor_df:", paste(dim(predictor_df), collapse = " x ")))
    
    # Run the Tweedieverse function
    write_log("Running Tweedieverse...")
    tryCatch({
      Tweedieverse::Tweedieverse(
        input_features = as.data.frame(count_table),
        input_metadata = predictor_df,
        output = output_dir,
        base_model = mymethod,
        abd_threshold = 0,
        prev_threshold = 0.0,
        var_threshold = 0,
        entropy_threshold = 0,
        random_effects = subject
      )
      write_log("Tweedieverse completed.")
    }, error = function(e) {
      write_log(paste("Error in Tweedieverse:", e$message))
      return(NULL)
    })
    
    # Check if the results file was created
    if (!file.exists(results_file)) {
      write_log("Results file does not exist")
      stop("Results file does not exist")
    }
    write_log("Results file exists.")
    
    # Read the results
    myTable <- read_tsv(results_file, show_col_types = FALSE)
    myTable <- as.data.frame(myTable)
    
    # Debug: Check myTable
    write_log("myTable:")
    write_log(paste(capture.output(head(myTable)), collapse = "\n"))
    write_log(paste("Dimensions of myTable:", paste(dim(myTable), collapse = " x ")))
    
    # Ensure myTable has the expected structure
    if (nrow(myTable) == 0 || !all(c("feature", "pval") %in% colnames(myTable))) {
      write_log("Results table has unexpected structure")
      stop("Results table has unexpected structure")
    }
    
    # Define a result data frame
    result_df <- data.frame(
      Feature = myTable$feature,
      pval = myTable$pval,
      pval.adj = p.adjust(myTable$pval, method = "fdr"),
      Method = "CPLM",
      stringsAsFactors = FALSE
    )
    rownames(result_df) <- result_df$Feature
    
    write_log("Returning result_df")
    write_log(paste(capture.output(head(result_df)), collapse = "\n"))
    
    return(result_df)
  }
  
  # Run testDA with cores set to 1
  final <- DAtest::testDA(
    data = df,
    predictor = vec,
    tests = c("zzz"),
    args = list(zzz = list(FUN = CPLM)),
    cores = 1,
    paired = subject)
  
  final$table = final$data
  
}

write_tsv(as.data.frame(final$table),opt$output)
