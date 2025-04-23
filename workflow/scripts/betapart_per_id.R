#!/usr/local/bin/Rscript --vanilla



set.seed(1234)
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(require(tidyverse)))
suppressWarnings(suppressMessages(library(devtools)))
suppressWarnings(suppressMessages(load_all(path = "workflow/scripts/betapart")))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="tsv file", metavar="Features input file formatted as tsv from a biom file"),
  make_option(c("-z", "--diffoutput"), type="character", default=NULL, help="Ingroup variation Output file", metavar="Diff Output file"),
  make_option(c("-o", "--permuoutput"), type="character", default=NULL, help="output folder to save the permutation files", metavar="Output permutation files"),
  make_option(c("-r", "--reps"), type="character", default=10, help="number of samples per each permutation", metavar="Samples per permutation"),
  make_option(c("-p", "--perm"), type="character", default=1000, help="number of permutations to run", metavar="Number of permutations"),
  make_option(c("-d", "--distance"), type="character", default=NULL, help="Type of distance matrix. Either the string 'bray' or 'jaccard'.", metavar="Distance matrix")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

biomFile = opt$input
output = opt$diffoutput 
distance = opt$distance
numOfSamplesPerTrial = opt$reps
NumOfPermutations = opt$perm
outputPerm = opt$permuoutput

subsetTable = suppressMessages(read.delim(biomFile, row.names = 1))


if (distance == "bray"){
  phylo.popu = beta.sample.abund(subsetTable, index.family = "bray", sites= as.numeric(numOfSamplesPerTrial), samples=as.numeric(NumOfPermutations))
  populationTable = phylo.popu$sampled.values
} else if (distance == "jaccard") {
  phylo.popu = beta.sample(vegan::decostand(x=subsetTable, method="pa"), index.family = "jaccard", sites= as.numeric(numOfSamplesPerTrial), samples=as.numeric(NumOfPermutations))  
  populationTable = phylo.popu$sampled.values
}


sink(output)

if (distance == "bray"){
print(paste("# Per Baselga (2017), Bray-Curtis dissimilarity can be broken down to 'Balanced variation in abundance' and 'Abundance gradients'. Simply put, 'balanced variation in abundance' is analogous to 'Spatial turnover' in incidence terms, that is, replacement by another species. 'abundance gradients' means that one sample is a subset of the other. This is analogous to 'nestness' in incidence terms, that is, removal and no replacement by another species (one sample is a subset of the other). Raw Bray-Curtis is a combination of these two differences. This script will calculate raw Bray-Curtis dissimilarity, as well as the two components of Bray-Curtis."))

brayCalc.core = betapart.core.abund(subsetTable)
brayCalc.multi = beta.multi.abund(brayCalc.core) # multiple site measures (this will give you the overall values)
print(paste0("Within group values: Bray Curtis: ",brayCalc.multi$beta.BRAY," Balanced: ", brayCalc.multi$beta.BRAY.BAL," Gradients: ", brayCalc.multi$beta.BRAY.GRA))
print(paste0("This means that Balanced represents:", (brayCalc.multi$beta.BRAY.BAL/brayCalc.multi$beta.BRAY)*100,"% and Gradients represents ", (brayCalc.multi$beta.BRAY.GRA/brayCalc.multi$beta.BRAY)*100, "%"))

} else if (distance == "jaccard"){
  
  print(paste("# Per Baselga (2012), Jaccard dissimilarity can be broken down to 'Turn-over' (replacement) and 'Nestedness'. Simply put, 'Turn-over' is replacement of one feature by another. 'Nestedness' means that one sample is a subset of the other. Raw Jaccard distances are a combination of these two types of differences. This script will calculate raw Jaccard dissimilarity, as well as the two components of Jaccard."))
  jacCalc.core = betapart.core(vegan::decostand(x=subsetTable, method="pa"))
  jacCalc.multi = beta.multi(jacCalc.core, index.family="jac") # multiple site measures (this will give you the overall values)
  print(paste0("Within group values: Jaccard: ", jacCalc.multi$beta.JAC," Turn-over: ", jacCalc.multi$beta.JTU," Nestedness: ", jacCalc.multi$beta.JNE))
  print(paste0("This means that Turn-over represents:", (jacCalc.multi$beta.JTU/jacCalc.multi$beta.JAC)*100,"% and Gradients represents ", (jacCalc.multi$beta.JNE/jacCalc.multi$beta.JAC)*100, "%"))
    
}

print(paste0("Permutation-based mean values for balanced, gradient, and full matrix of ", distance))
print(paste(phylo.popu$mean.values))
print(paste("Permutation-based SD values for balanced, gradient, and full matrix of", distance))
print(paste(phylo.popu$sd.values))

sink()

write.table(populationTable,outputPerm,sep="\t",row.names = FALSE)