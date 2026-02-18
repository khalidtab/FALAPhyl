#!/usr/local/bin/Rscript --vanilla
#
# plot_dunns_heatmap.R
#
# Creates a heatmap of Dunn's Z-scores across multiple beta diversity
# distance metrics for a given grouping variable. Each cell shows the
# Z-score for a pairwise group comparison, colored on a diverging scale.
# Significance annotations (* p<0.05, ** p<0.01, *** p<0.001) are
# overlaid on each cell.
#
# Input: Directory containing Dunn's test output files
#        (Dunns–{dist}–{group}.txt from the dunns_ rule)
# Output: SVG heatmap figure
#
# Usage:
#   Rscript --vanilla plot_dunns_heatmap.R \
#       -i data/distance/PerGroupDistCompare/{sample}/ \
#       -g {group} -o output.svg -x 12 -y 8

suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(svglite)))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input directory containing Dunn's test result files"),
  make_option(c("-g", "--group"), type = "character", default = NULL,
              help = "Grouping variable name (to match filenames)"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output SVG file path"),
  make_option(c("-x", "--width"), type = "numeric", default = 12,
              help = "Plot width in inches [default %default]"),
  make_option(c("-y", "--height"), type = "numeric", default = 8,
              help = "Plot height in inches [default %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$group) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Required arguments: -i (input dir), -g (group), -o (output file)")
}

input_dir <- opt$input
group_var <- opt$group
output_file <- opt$output

# Find all Dunn's result files for this group
pattern <- paste0("Dunns.+", group_var, "\\.txt$")
dunn_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)

if (length(dunn_files) == 0) {
  stop("No Dunn's test files found matching pattern '", pattern, "' in ", input_dir)
}

message("Found ", length(dunn_files), " Dunn's test result files.")

# Parse each file and extract the tabular Dunn's results
all_results <- list()

for (f in dunn_files) {
  # Extract distance metric from filename: Dunns–{dist}–{group}.txt
  fname <- basename(f)
  # Handle en-dash (–) or hyphen (-) separators
  dist_metric <- sub("^Dunns[–-](.+)[–-][^–-]+\\.txt$", "\\1", fname)

  # Read the file — it may have header lines from sink() before the table
  lines <- readLines(f, warn = FALSE)

  # Find the header line of the Dunn's table (starts with "Category.1")
  header_idx <- grep("^Category\\.1", lines)

  if (length(header_idx) == 0) {
    message("Skipping ", fname, " — no Dunn's table found (may be 2-group or non-significant).")
    next
  }

  # Read from the header line onwards
  df <- tryCatch(
    read.delim(textConnection(lines[header_idx[1]:length(lines)]),
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      message("Skipping ", fname, " — parse error: ", e$message)
      return(NULL)
    }
  )

  if (is.null(df) || nrow(df) == 0) next

  df$distance <- dist_metric
  all_results[[fname]] <- df
}

if (length(all_results) == 0) {
  stop("No valid Dunn's test results found. All files may be 2-group comparisons or non-significant.")
}

results <- bind_rows(all_results)

# Create comparison label
results$comparison <- paste0(results$Category.1, " vs. ", results$Category.2)

# Significance stars
results$sig <- case_when(
  results$P.adjusted < 0.001 ~ "***",
  results$P.adjusted < 0.01  ~ "**",
  results$P.adjusted < 0.05  ~ "*",
  TRUE                       ~ ""
)

# Z-score label (rounded, with significance)
results$label <- paste0(round(results$Z, 2), results$sig)

# Determine symmetric Z-score limits for diverging color scale
z_limit <- max(abs(results$Z), na.rm = TRUE) * 1.05

# Order comparisons and distances for consistent layout
results$comparison <- factor(results$comparison,
                             levels = sort(unique(results$comparison)))
results$distance <- factor(results$distance,
                           levels = sort(unique(as.character(results$distance))))

# Build heatmap
p <- ggplot(results, aes(x = distance, y = comparison, fill = Z)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.2, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-z_limit, z_limit),
    name = "Dunn's Z"
  ) +
  labs(
    x = "Distance metric",
    y = "Pairwise comparison",
    caption = "* p < 0.05, ** p < 0.01, *** p < 0.001 (BH-adjusted)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    plot.caption = element_text(size = 9, hjust = 0),
    legend.position = "right"
  )

# Save
ggsave(output_file, plot = p, width = opt$width, height = opt$height, device = svglite)
message("Heatmap saved to: ", output_file)
