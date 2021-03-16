suppressMessages(library(tidyverse))
suppressMessages(library("optparse"))
suppressMessages(library(reshape2))


option_list = list(
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="The mapping file", metavar="Mapping file"),
  make_option(c("-i", "--input"), type="character", default=NULL, help="The alpha diversity file", metavar="Alpha diversity file"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="The groups column in your mapping file you are comparing", metavar="Groups"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Path to output file of the non-parametric test", metavar="Path to output file")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

# Load files
map = opt$mapping
alpha = opt$input
group = opt$group
output = opt$output

map = suppressMessages(read_tsv(map))
alpha = suppressMessages(read_tsv(alpha))

# Merge then keep only the pertinent columns
myNames = c("#SampleID","alpha")
colnames(alpha) = myNames
merged = merge(map,alpha) %>% .[, c(myNames, group)]
colnames(merged)[3] = "group"


#merged = filter(merged, group == c("normal","overweight"))

if (length(unique(merged$group)) > 2){
  test = kruskal.test(merged$alpha ~ merged$group, data = merged)
  posthoc = pairwise.wilcox.test(merged$alpha, merged$group, p.adjust.method = "fdr")
  melted = melt(posthoc$p.value, na.rm = TRUE)
  method = paste0(posthoc$method,". Post-hoc testing method: ",posthoc$p.adjust.method)
  sink(output)
  print(test)
  print(method)
  print(melted)
  sink()
  } else {
  test = wilcox.test(merged$alpha ~ merged$group, data = merged)
  sink(output)
  print(test)
  sink()
}



