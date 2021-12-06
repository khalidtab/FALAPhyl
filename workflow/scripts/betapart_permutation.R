#!/usr/local/bin/Rscript --vanilla



set.seed(1234)

suppressWarnings(suppressMessages(library(devtools)))	
suppressWarnings(suppressMessages(load_all(path = "workflow/scripts/betapart")))
suppressWarnings(suppressMessages(library(optparse)))	
suppressWarnings(suppressMessages(require(tidyverse)))	
suppressWarnings(suppressMessages(require(vegan)))	
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(phyloseq)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="tsv file", metavar="Features input file formatted as tsv from a biom file"),
  make_option(c("-r", "--reps"), type="character", default=10, help="number of samples per each permutation", metavar="Samples per permutation"),
  make_option(c("-p", "--perm"), type="character", default=1000, help="number of permutations to run", metavar="Number of permutations"),
  make_option(c("-d", "--distance"), type="character", default=NULL, help="Type of distance matrix. Either the string 'bray' or 'jaccard'.", metavar="Distance matrix"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder to save the files", metavar="Output folder")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

subsetTable = opt$input
output = opt$output
numOfSamplesPerTrial = opt$reps
NumOfPermutations = opt$perm
distance = opt$distance

myBasename = basename(subsetTable)
subsetTable = suppressMessages(read.delim(subsetTable, row.names = 1)) %>% t(.)

if (distance == "bray"){
phylo.popu = beta.sample.abund(subsetTable, index.family = "bray", sites= as.numeric(numOfSamplesPerTrial), samples=as.numeric(NumOfPermutations))
populationTable = phylo.popu$sampled.values
} else if (distance == "jaccard") {
  phylo.popu = beta.sample(vegan::decostand(x=subsetTable, method="pa"), index.family = "jaccard", sites= as.numeric(numOfSamplesPerTrial), samples=as.numeric(NumOfPermutations))  
  populationTable = phylo.popu$sampled.values
}
sink(paste0(output,"/betapart_permu_",myBasename,"_",distance,"_mean_sd.txt"))

print(paste0(myBasename,": Permutation-based mean values for balanced, gradient, and full matrix of ", distance))
print(paste(phylo.popu$mean.values))
print(paste(myBasename,": Permutation-based SD values for balanced, gradient, and full matrix of", distance))
print(paste(phylo.popu$sd.values))


sink()

write.table(populationTable,paste0(output,"permutations/betapart_permutations_",distance,"_",myBasename),sep="\t",row.names = FALSE)


