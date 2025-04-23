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

all_comparisons <- unique(myFile$comparison)
results_list <- list()

for (comp in all_comparisons) {
  sub_data <- myFile %>% filter(comparison == comp)
  n_groups <- length(unique(sub_data$test_type))
  
  if (n_groups < 2) {
    message(paste("Skipping:", comp, "— only one group"))
    results_list[[comp]] <- data.frame(Comparison = comp, Message = "Only one group")
    
  } else if (n_groups == 2) {
    # Wilcoxon test
    test_result <- wilcox.test(score ~ test_type, data = sub_data)
    results_list[[comp]] <- data.frame(
      Comparison = comp,
      Method = "Wilcoxon",
      p_value = test_result$p.value
    )
    
  } else {
    # Dunn test
    dt <- dunn.test(sub_data$score, sub_data$test_type, method = "bh", kw = FALSE)
    df <- data.frame(
      Comparison = comp,
      ComparisonGroup = dt$comparisons,
      Z = dt$Z,
      P = dt$P,
      P.adjusted = dt$P.adjusted
    )
    results_list[[comp]] <- df
  }
}

# Combine all into a single data frame
all_results <- bind_rows(results_list)

output <- capture.output(print(all_results))
writeLines(output, opt$output)
