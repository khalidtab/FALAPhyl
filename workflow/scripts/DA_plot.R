#!/usr/local/bin/Rscript --vanilla

library("optparse")
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
library("grid")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Folder with the files to create the graph", metavar="Input folder"),
  make_option(c("-a", "--auc"), type="character", default=NULL, help="Path for AUC graph", metavar="Path for AUC graph"),
  make_option(c("-f", "--fdr"), type="character", default=NULL, help="Path for FDR graph", metavar="Path for FDR graph"),
  make_option(c("-p", "--power"), type="character", default=NULL, help="Path for Power graph", metavar="Path for Power graph"),
  make_option(c("-s", "--score"), type="character", default=NULL, help="Path for Score graph", metavar="Path for Scores graph"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output tsv file", metavar="Output tsv File name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

print(opt$input)

# List and flatten file paths
myFiles = list.dirs(path = opt$input) %>%
  grep("effectSize", x = ., value = TRUE) %>%
  sapply(function(x) list.files(x, full.names = TRUE)) %>%
  unlist()

# Extract comparisons (e.g., "Subgingival_plaque!Supragingival_plaque")
comparisons = basename(myFiles) %>%
  gsub(pattern = "\\.tsv$", replacement = "", .) %>%
  unique()

# Extract types of test (e.g., "abc", "wil")
typesOfTests = gsub(pattern = ".*–[0-9]*", replacement = "", x = myFiles) %>%
  gsub(pattern = "/.*", replacement = "", .) %>%
  gsub(pattern = "[0-9]–", replacement = "", .) %>%
  unique()

# Extract effect sizes from path
get_effect_size <- function(path) {
  stringr::str_extract(path, "effectSize–\\d+") %>%
    gsub("effectSize–", "", .) %>%
    as.numeric()
}

# Extract comparisons from file names
get_comparison <- function(path) {
  tools::file_path_sans_ext(basename(path))
}

# Function to read all tables for a specific test
theTestTables = function(theTest){
  files = grep(theTest, myFiles, value = TRUE)
  bind_rows(lapply(files, function(x) {
    read_tsv(x, show_col_types = FALSE) %>%
      mutate(
        comparison = get_comparison(x),
        effect_size = get_effect_size(x),
        test_type = theTest
      )
  }))
}

# Build the final table
theTestFiles = lapply(typesOfTests, theTestTables) %>%
  bind_rows() %>%
  mutate(score = (AUC - 0.5) * Power - FDR)


###### Functions ######
plotAUC = function(){

# Ensure effect size is factor for legend/transparency
theTestFiles <- theTestFiles %>%
  mutate(effect_size_factor = as.factor(effect_size))

# Create numeric Method for x-axis manipulation
method_levels <- unique(theTestFiles$Method)
method_map <- data.frame(Method = method_levels, Method_numeric = 1:length(method_levels))

# Join numeric mapping back into main data
theTestFiles <- theTestFiles %>%
  left_join(method_map, by = "Method")

# Global means per method, to be nudged right
global_meansAUC <- theTestFiles %>%
  group_by(Method, Method_numeric) %>%
  summarise(GlobalMean = mean(AUC, na.rm = TRUE), .groups = "drop") %>%
  mutate(Method_nudged = Method_numeric + 0.4)

# Effect-size-specific means per method & comparison, to be nudged + jittered
effect_meansAUC <- theTestFiles %>%
  group_by(comparison, Method, Method_numeric, effect_size_factor) %>%
  summarise(EffectMean = mean(AUC, na.rm = TRUE), .groups = "drop") %>%
  mutate(Method_nudged = Method_numeric + 0.4)

# Plot
AUC <- ggplot(theTestFiles, aes(x = Method_numeric, y = AUC)) +
  
  # Raw AUC points
  geom_jitter(aes(alpha = effect_size_factor), color = "black", width = 0.2, size = 1) +
  
  # Effect size means: white circles, nudged + jittered horizontally
  geom_point(data = effect_meansAUC,
             aes(x = Method_nudged, y = EffectMean, alpha = effect_size_factor),
             inherit.aes = FALSE,
             shape = 21, fill = "blue", color = "black", size = 5,
             position = position_jitter(width = 0.1, height = 0)) +
  
  # Global means: black triangle, nudged right (no jitter)
  geom_point(data = global_meansAUC,
             aes(x = Method_nudged, y = GlobalMean),
             inherit.aes = FALSE,
             shape = 24, fill = "black", color = "white", size = 3) +
  
  # Reference line at AUC = 0.5
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  
  # Facet per comparison
  facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
  
  # Use numeric scale for x, with original method labels
  scale_x_continuous(
    breaks = method_map$Method_numeric,
    labels = method_map$Method
  ) +
  
  # Labels & styling
  labs(
    x = "Method",
    y = "AUC",
    alpha = "Effect Size",
    caption = paste(
      "• Black dots: individual AUCs (transparency = effect size)",
      "• Blue circles: per-effect-size means (transparency = effect size)",
      "• Black triangles: global mean across comparisons (nudged)",
      "• Dashed red line = AUC = 0.5 = random performance, means spiked features are randomly spread with non-spiked.",
      "Therefore, we want an AUC as high as possible. That is, spiked features should have low p-values (ie they have been identified)",
      sep = "\n"
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold", size = 11),
    plot.caption = element_text(size = 9),
    legend.position = "right"
  )
return(AUC)
}
plotPower = function(){
  #### POWER ####
  # Ensure effect size is factor for alpha/transparency
  theTestFiles <- theTestFiles %>%
    mutate(effect_size_factor = as.factor(effect_size))
  
  # Create numeric mapping for Method
  method_levels <- unique(theTestFiles$Method)
  method_map <- data.frame(Method = method_levels, Method_numeric = 1:length(method_levels))
  
  # Join numeric x-position to data
  theTestFiles <- theTestFiles %>%
    left_join(method_map, by = "Method")
  
  # Global mean Power per method
  global_meansPower <- theTestFiles %>%
    group_by(Method, Method_numeric) %>%
    summarise(GlobalMean = mean(Power, na.rm = TRUE), .groups = "drop") %>%
    mutate(Method_nudged = Method_numeric + 0.4)
  
  # Effect-size-specific mean Power
  effect_meansPower <- theTestFiles %>%
    group_by(comparison, Method, Method_numeric, effect_size_factor) %>%
    summarise(EffectMean = mean(Power, na.rm = TRUE), .groups = "drop") %>%
    mutate(Method_nudged = Method_numeric + 0.4)
  
  # Plot
  Power <- ggplot(theTestFiles, aes(x = Method_numeric, y = Power)) +
    
    # Raw Power points
    geom_jitter(aes(alpha = effect_size_factor), color = "black", width = 0.2, size = 1) +
    
    # Effect size means: green circles, jittered horizontally
    geom_point(data = effect_meansPower,
               aes(x = Method_nudged, y = EffectMean, alpha = effect_size_factor),
               inherit.aes = FALSE,
               shape = 21, fill = "darkgreen", color = "black", size = 5,
               position = position_jitter(width = 0.1, height = 0)) +
    
    # Global means: black triangle, nudged only
    geom_point(data = global_meansPower,
               aes(x = Method_nudged, y = GlobalMean),
               inherit.aes = FALSE,
               shape = 24, fill = "black", color = "white", size = 3) +
    
    # Facet by comparison
    facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
    
    # Numeric x-axis with original method names
    scale_x_continuous(
      breaks = method_map$Method_numeric,
      labels = method_map$Method
    ) +
    
    # Labels & caption
    labs(
      x = "Method",
      y = "Power",
      alpha = "Effect Size",
      caption = paste(
        "• Black dots: individual Power values (transparency = effect size)",
        "• Green circles: per-effect-size means (jittered + nudged)",
        "• Black triangles: global mean across all comparisons (nudged)",
        "• Power is the proportion of spiked features that are significant after multiple p-value corrections.",
        "The higher the power, the better.",
        sep = "\n"
      )
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 9),
      legend.position = "right"
    )
  
  return(Power)
}
plotFDR = function(){
  #### FDR ####
  # Ensure effect size is a factor for alpha/transparency
  theTestFiles <- theTestFiles %>%
    mutate(effect_size_factor = as.factor(effect_size))
  
  # Numeric x-axis for method
  method_levels <- unique(theTestFiles$Method)
  method_map <- data.frame(Method = method_levels, Method_numeric = 1:length(method_levels))
  
  theTestFiles <- theTestFiles %>%
    left_join(method_map, by = "Method")
  
  # Global FDR means per method
  global_meansFDR <- theTestFiles %>%
    group_by(Method, Method_numeric) %>%
    summarise(GlobalMean = mean(FDR, na.rm = TRUE), .groups = "drop") %>%
    mutate(Method_nudged = Method_numeric + 0.4)
  
  # Effect-size-specific FDR means
  effect_meansFDR <- theTestFiles %>%
    group_by(comparison, Method, Method_numeric, effect_size_factor) %>%
    summarise(EffectMean = mean(FDR, na.rm = TRUE), .groups = "drop") %>%
    mutate(Method_nudged = Method_numeric + 0.4)
  
  # Plot
  FDR <- ggplot(theTestFiles, aes(x = Method_numeric, y = FDR)) +
    
    # Jittered individual FDR points
    geom_jitter(aes(alpha = effect_size_factor), color = "black", width = 0.2, size = 1) +
    
    # Effect size means: white circles, jittered horizontally
    geom_point(data = effect_meansFDR,
               aes(x = Method_nudged, y = EffectMean, alpha = effect_size_factor),
               inherit.aes = FALSE,
               shape = 21, fill = "maroon", color = "black", size = 5,
               position = position_jitter(width = 0.1, height = 0)) +
    
    # Global means: black triangle, nudged right
    geom_point(data = global_meansFDR,
               aes(x = Method_nudged, y = GlobalMean),
               inherit.aes = FALSE,
               shape = 24, fill = "black", color = "white", size = 3) +
    
    # Facet per comparison
    facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
    
    # Method names on x-axis
    scale_x_continuous(
      breaks = method_map$Method_numeric,
      labels = method_map$Method
    ) +
    
    # Labels and caption
    labs(
      x = "Method",
      y = "False Discovery Rate (FDR)",
      alpha = "Effect Size",
      caption = paste(
        "• Black dots: individual FDR values (transparency = effect size)",
        "• Maroon circles: per-effect-size means (jittered + nudged)",
        "• Black triangles: global mean across comparisons (nudged)",
        "• FDR is the proportion of significant features that were NOT spiked (after multiple comparison adjustment) that were not spiked and therefore shouldn't be significant",
        "• Lower is better — high FDR indicates false positives",
        sep = "\n"
      )
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 9),
      legend.position = "right"
    )
  
  return(FDR)
}
plotScore = function(){
  #### SCORE ####
  # Ensure effect size is a factor for alpha/transparency
  theTestFiles <- theTestFiles %>%
    mutate(effect_size_factor = as.factor(effect_size))
  
  # Numeric mapping for Method axis
  method_levels <- unique(theTestFiles$Method)
  method_map <- data.frame(Method = method_levels, Method_numeric = 1:length(method_levels))
  
  theTestFiles <- theTestFiles %>%
    left_join(method_map, by = "Method")
  
  # Global mean scores per method
  global_meansScore <- theTestFiles %>%
    group_by(Method, Method_numeric) %>%
    summarise(GlobalMean = mean(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(Method_nudged = Method_numeric + 0.4)
  
  # Per-effect-size mean scores
  effect_meansScore <- theTestFiles %>%
    group_by(comparison, Method, Method_numeric, effect_size_factor) %>%
    summarise(EffectMean = mean(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(Method_nudged = Method_numeric + 0.4)
  
  # Plot
  Score <- ggplot(theTestFiles, aes(x = Method_numeric, y = score)) +
    
    # Raw score points
    geom_jitter(aes(alpha = effect_size_factor), color = "black", width = 0.2, size = 1) +
    
    # Effect size means (blue circles, jittered horizontally)
    geom_point(data = effect_meansScore,
               aes(x = Method_nudged, y = EffectMean, alpha = effect_size_factor),
               inherit.aes = FALSE,
               shape = 21, fill = "orange", color = "black", size = 5,
               position = position_jitter(width = 0.1, height = 0)) +
    
    # Global mean scores (black triangles, nudged right)
    geom_point(data = global_meansScore,
               aes(x = Method_nudged, y = GlobalMean),
               inherit.aes = FALSE,
               shape = 24, fill = "black", color = "white", size = 3) +
    
    # Facet by comparison
    facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
    
    # Original Method labels on x-axis
    scale_x_continuous(
      breaks = method_map$Method_numeric,
      labels = method_map$Method
    ) +
    
    # Labels & styling
    labs(
      x = "Method",
      y = "Score",
      title = "Score",
      alpha = "Effect Size",
      caption = paste(
        "• Black dots: individual scores (transparency = effect size)",
        "• Orange circles: per-effect-size mean scores (jittered + nudged)",
        "• Black triangles: global mean score across comparisons (nudged)",
        "• Score = (AUC - 0.5) * Power - FDR",
        "• The higher the Score, the better the method is estimated to be",
        sep = "\n"
      )
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 9),
      legend.position = "right"
    )
  
  return(Score)
}
########################


AUC = plotAUC()
Power = plotPower()
FDR = plotFDR()
Score = plotScore()

rm(plotAUC,plotPower,plotFDR,plotScore)
gc()

write_tsv(theTestFiles,opt$output)

svg(opt$auc, width = 10, height = (10*length(comparisons)))
AUC
dev.off()

svg(opt$power,width = 10, height = (10*length(comparisons)))
Power
dev.off()

svg(opt$fdr,width = 10, height = (10*length(comparisons)))
FDR
dev.off()

svg(opt$score ,width = 10, height = (10*length(comparisons)))
Score
dev.off()