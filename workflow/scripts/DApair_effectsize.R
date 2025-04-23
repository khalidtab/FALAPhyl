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
  make_option(c("-p", "--tmp"), type="character", default=NULL, help="Path to temporary folder", metavar="Path to temporary folder"),
  make_option(c("-u", "--subjectID"), type="character", default=NULL, help="Column name of the subject ID in the mapping file", metavar="Column name of the subject ID in the mapping file"),  
  make_option(c("-l", "--log"), type="character", default=NULL, help="Log file name", metavar="Log file name"),
  make_option(c("-e", "--effectsize"), type="character", default=NULL, help="effect size", metavar="effect size")
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
df = t(df) %>% as.data.frame(.)

# Load mapping file
map = opt$mapping
map = read.csv(map,sep="\t") %>% as.data.frame(.)

category = opt$category
catNum = which(colnames(map) == category)
working_map = cbind(as.character(map[,1]),
                    as.character(map[,catNum])) %>% as.data.frame(.)
colnames(working_map) = c("SampleID","condition")
vec = working_map$condition %>% as.factor(.)

log_file = opt$log

tmp_folder = opt$tmp

mymethod = opt$thetest

effectSize = as.numeric(opt$effectsize)

subjectID = opt$subjectID

if (mymethod == "abc"){
  final=testDA(df, predictor = vec,tests=c("kru","abc"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "adx") {
  final=testDA(df, predictor = vec,tests=c("kru","adx"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "aov") {
  final=testDA(df, predictor = vec,tests=c("kru","aov"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "lao") {
  final=testDA(df, predictor = vec,tests=c("kru","lao"), effectSize = effectSize, cores = 1)
  
  
}  else if (mymethod == "lao2") {
  final=testDA(df, predictor = vec,tests=c("kru","lao2"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "bay") {
  final=testDA(df, predictor = vec,tests=c("kru","bay"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "pea") {
  final=testDA(df, predictor = vec,tests=c("kru","pea"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "spe") {
  final=testDA(df, predictor = vec,tests=c("kru","spe"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "ds2x") {
  final=testDA(df, predictor = vec,tests=c("kru","ds2x"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "ds2") {
  final=testDA(df, predictor = vec,tests=c("kru","ds2"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "ere") {
  final=testDA(df, predictor = vec,tests=c("kru","ere"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "ere2") {
  final=testDA(df, predictor = vec,tests=c("kru","ere2"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "erq") {
  final=testDA(df, predictor = vec,tests=c("kru","erq"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "erq2") {
  final=testDA(df, predictor = vec,tests=c("kru","erq2"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "fri") {
  final=testDA(df, predictor = vec,tests=c("kru","fri"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "neb") {
  final=testDA(df, predictor = vec,tests=c("kru","neb"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "poi") {
  final=testDA(df, predictor = vec,tests=c("kru","poi"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "qpo") {
  final=testDA(df, predictor = vec,tests=c("kru","qpo"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "znb") {
  final=testDA(df, predictor = vec,tests=c("kru","znb"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "zpo") {
  final=testDA(df, predictor = vec,tests=c("kru","zpo"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "lim") {
  final=testDA(df, predictor = vec,tests=c("kru","lim"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "lli") {
  final=testDA(df, predictor = vec,tests=c("kru","lli"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "lli2") {
  final=testDA(df, predictor = vec,tests=c("kru","lli2"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "vli") {
  final=testDA(df, predictor = vec,tests=c("kru","vli"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "lrm") {
  final=testDA(df, predictor = vec,tests=c("kru","lrm"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "llm") {
  final=testDA(df, predictor = vec,tests=c("kru","llm"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "llm2") {
  final=testDA(df, predictor = vec,tests=c("kru","llm2"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "msf") {
  final=testDA(df, predictor = vec,tests=c("kru","msf"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "zig") {
  final=testDA(df, predictor = vec,tests=c("kru","zig"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "mva") {
  final=testDA(df, predictor = vec,tests=c("kru","mva"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "per") {
  final=testDA(df, predictor = vec,tests=c("kru","per"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "qua") {
  final=testDA(df, predictor = vec,tests=c("kru","qua"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "sam") {
  final=testDA(df, predictor = vec,tests=c("kru","sam"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "ttt") {
  final=testDA(df, predictor = vec,tests=c("kru","ttt"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "ltt") {
  final=testDA(df, predictor = vec,tests=c("kru","ltt"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "ltt2") {
  final=testDA(df, predictor = vec,tests=c("kru","ltt2"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "ttr") {
  final=testDA(df, predictor = vec,tests=c("kru","ttr"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "wil") {
  final=testDA(df, predictor = vec,tests=c("kru","wil"), effectSize = effectSize, cores = 1)
  
  
}  else if (mymethod == "ttc") {
  final=testDA(df, predictor = vec,tests=c("kru","ttc"), effectSize = effectSize, cores = 1)
  
  
}  else if (mymethod == "lic") {
  final=testDA(df, predictor = vec,tests=c("kru","lic"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "tta") {
  final=testDA(df, predictor = vec,tests=c("kru","tta"), effectSize = effectSize, cores = 1)
  
  
} else if (mymethod == "lia") {
  final=testDA(df, predictor = vec,tests=c("kru","lia"), effectSize = effectSize, cores = 1)
  
} else if (mymethod %in% c("CPLM","ZICP","ZSCP","ZACP")){

  library("Tweedieverse")
  
 df <- apply(df, 2, function(x) x / sum(x))
  
  CPLM <- function(count_table, predictor, paired, covars) {
    
    # Create a unique output directory based on timestamp
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- paste0(tmp_folder, timestamp, "_CPLM")
    dir.create(output_dir, recursive = TRUE)
    results_file <- paste0(output_dir, "/all_results.tsv")
    
    print("Entering function")
    
    # Debug: Print inputs
    print("count_table:")
    print(paste(capture.output(head(count_table)), collapse = "\n"))
    print(paste("Dimensions of count_table:", paste(dim(count_table), collapse = " x ")))
    
    print("predictor:")
    print(paste(capture.output(head(predictor)), collapse = "\n"))
    print(paste("Length of predictor:", length(predictor)))
    
    # Prepare the predictor data frame
    predictor_df <- data.frame(sampleid = colnames(count_table), metadata = as.character(predictor))
    rownames(predictor_df) <- predictor_df$sampleid
    predictor_df$sampleid <- NULL
    
    print("predictor_df:")
    print(paste(capture.output(head(predictor_df)), collapse = "\n"))
    print(paste("Dimensions of predictor_df:", paste(dim(predictor_df), collapse = " x ")))
    
    # Run the Tweedieverse function
    print("Running Tweedieverse...")
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
      )
      print("Tweedieverse completed.")
    }, error = function(e) {
      print(paste("Error in Tweedieverse:", e$message))
      return(NULL)
    })
    
    # Check if the results file was created
    if (!file.exists(results_file)) {
      print("Results file does not exist")
      stop("Results file does not exist")
    }
    print("Results file exists.")
    
    # Read the results
    myTable <- read_tsv(results_file, show_col_types = FALSE)
    myTable <- as.data.frame(myTable)
    
    # Debug: Check myTable
    print("myTable:")
    print(paste(capture.output(head(myTable)), collapse = "\n"))
    print(paste("Dimensions of myTable:", paste(dim(myTable), collapse = " x ")))
    
    # Ensure myTable has the expected structure
    if (nrow(myTable) == 0 || !all(c("feature", "pval") %in% colnames(myTable))) {
      print("Results table has unexpected structure")
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
    
    print("Returning result_df")
    print(paste(capture.output(head(result_df)), collapse = "\n"))
    
    return(result_df)
  }
  
  # Run testDA with cores set to 1
  final_results <- DAtest::testDA(
    data = df,
    predictor = vec,
    tests = c("zzz"),
    args = list(zzz = list(FUN = CPLM)), relative = FALSE, cores = 1, effectSize = effectSize)

  library(ROCR)
  
  # Assuming 'data' is the input data frame
  data <- as.data.frame(final_results$results)
  
  # Identify the columns for adjusted p-values and spiked status
  pval_cols <- grep("zzz.pval.adj", colnames(data), value = TRUE)
  spiked_cols <- grep("zzz.Spiked", colnames(data), value = TRUE)
  raw_pval_cols <- gsub("zzz.pval.adj", "zzz.pval", pval_cols)
  
  # Define the significance threshold
  significance_threshold <- 0.05
  
  # Function to calculate metrics for each run
  calculate_metrics <- function(i) {
    # Extract p-values and spiked status for the current run
    pvals <- data[[pval_cols[i]]]
    raw_pvals <- data[[raw_pval_cols[i]]]
    spiked <- ifelse(data[[spiked_cols[i]]] == "Yes", 1, 0)
    
    # Calculate AUC
    pred <- prediction(-pvals, spiked)
    perf <- performance(pred, "auc")
    auc_value <- as.numeric(perf@y.values)
    
    # Calculate FPR
    non_spiked <- spiked == 0
    significant_non_spiked <- sum(raw_pvals[non_spiked] < significance_threshold)
    total_non_spiked <- sum(non_spiked)
    fpr <- significant_non_spiked / total_non_spiked
    
    # Calculate FDR
    significant_features <- pvals < significance_threshold
    false_discoveries <- sum(significant_features & non_spiked)
    total_significant <- sum(significant_features)
    fdr <- if (total_significant == 0) 0 else false_discoveries / total_significant
    
    # Calculate Power
    spiked_features <- spiked == 1
    significant_spiked <- sum(pvals[spiked_features] < significance_threshold)
    total_spiked <- sum(spiked_features)
    power <- significant_spiked / total_spiked
    
    # Return results as a list
    return(c(AUC = auc_value, FPR = fpr, FDR = fdr, Power = power))
  }
  
  # Use lapply to calculate metrics for all runs and convert the result to a data frame
  metrics <- do.call(rbind, lapply(seq_along(pval_cols), calculate_metrics))
  
  # Create the final summary data frame
  final_summary <- data.frame(
    Method = mymethod,  # Replace with your actual method name if different
    metrics,
    Run = seq_len(nrow(metrics)))
  
  final = NULL
  final$table = final_summary
}


write_tsv(as.data.frame(final$table),opt$output)
