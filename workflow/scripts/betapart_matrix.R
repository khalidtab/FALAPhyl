#!/usr/local/bin/Rscript --vanilla

set.seed(1234)
suppressWarnings(suppressMessages(library(devtools)))
suppressWarnings(suppressMessages(load_all(path = "workflow/scripts/betapart")))



#install.packages("workflow/scripts/betapart_1.5.4.tar.gz", repos = NULL, type ="source")
#library("betapart")
suppressWarnings(suppressMessages(library(optparse)))	
suppressWarnings(suppressMessages(require(tidyverse)))	
suppressWarnings(suppressMessages(require(vegan)))	
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(phyloseq)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Biom file", metavar="Features input file formatted as biom"),
  make_option(c("-d", "--distance"), type="character", default=NULL, help="Type of distance matrix. Either the string 'bray' or 'jaccard'.", metavar="Distance matrix"),
  make_option(c("-f", "--full"), type="character", default=NULL, help="output file to save the Bray-Curtis or Jaccard tsv matrix", metavar="Output matrix"),
  make_option(c("-r", "--replace"), type="character", default=NULL, help="output file to save the turnover/balanced (ie replacement) tsv matrix", metavar="Output replacement"),
  make_option(c("-n", "--noreplace"), type="character", default=NULL, help="output file to save the nestedness/gradient (ie no replacement) tsv matrix", metavar="Output no replacement")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


biomFile = opt$input
replace = opt$replace
noreplace = opt$noreplace
full= opt$full
distance = opt$distance

####### Functions #######
breakdown = function(biomMatrix,distance){
  biomMatrix = biomMatrix %>% t() %>% as.matrix()
  
  if (distance == "bray"){
    breakdown = beta.pair.abund(biomMatrix, index.family="bray")
  
  replacement   = breakdown$beta.bray.bal %>% as.matrix() %>% as.data.frame() 
  noReplacement = breakdown$beta.bray.gra %>% as.matrix %>% as.data.frame() 
  fullDistance  = breakdown$beta.bray %>% as.matrix %>% as.data.frame() 
  
  } else if (distance == "jaccard"){
    breakdown = beta.pair(vegan::decostand(x=biomMatrix, method="pa"), index.family="jaccard") #decostand changes abundance to presence/absence
  
    replacement   = breakdown$beta.jtu %>% as.matrix() %>% as.data.frame() 
    noReplacement = breakdown$beta.jne %>% as.matrix %>% as.data.frame() 
    fullDistance  = breakdown$beta.jac %>% as.matrix %>% as.data.frame() 
    
    }
  
  myDists = list(fullDistance,replacement,noReplacement)
  
  return(myDists)
}
#########################

# Make the phyloseq merged file
myTable = suppressMessages(readr::read_tsv(biomFile)) %>% as.data.frame(.)
rownames(myTable) = myTable[,1]
myTable[,1] = NULL

# Write the full dissimilarity, replacement and no replacement matricesbiom
myDists = breakdown(myTable,distance) 

write.table(x=myDists[[1]],file=paste0(full), quote = FALSE, sep = "\t",col.names = NA)
write.table(x=myDists[[2]],file=paste0(replace), quote = FALSE, sep = "\t",col.names = NA)
write.table(x=myDists[[3]],file=paste0(noreplace), quote = FALSE, sep = "\t",col.names = NA)