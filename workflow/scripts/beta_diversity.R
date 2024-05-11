#!/usr/local/bin/Rscript --vanilla

set.seed(1234)
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(require(tidyverse)))
suppressWarnings(suppressMessages(require(vegan)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(phyloseq)))
suppressWarnings(suppressMessages(library(ape)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="tsv file", metavar="Features input file formatted as tsv from a biom file"),
  make_option(c("-d", "--dissimilarity"), type="character", default=NULL, help="Dissimilarity type", metavar="Wanted dissimilarity type based on the input file options"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder to save the files", metavar="Output folder")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

tsvfile = opt$input 
dissimilarity = opt$dissimilarity
output = opt$output 

# Format the input file
mytsv = read.table(tsvfile,sep="\t",comment.char = "@",header = TRUE,row.names = 1)
mytsv = mytsv %>% mutate_if(is.character,as.numeric)


if (dissimilarity == "PhILR"){
# Generate tree
  mytsv =  otu_table(mytsv, taxa_are_rows = TRUE)
  mytree = suppressWarnings(dist(mytsv,method="euclidean")) %>% hclust(.,method="ward.D") %>% as.phylo(.) 

phylo = suppressMessages(merge_phyloseq(mytsv,mytree))

# Test tree
isRooted = suppressWarnings(ape::is.rooted(phy_tree(phylo)))
isBinary = suppressWarnings(ape::is.binary(phy_tree(phylo)))
if(!isRooted){stop("Tree needs to be rooted.",call.=FALSE)} #if not rooted
if(!isBinary){phy_tree(phylo) = multi2di(phylo@phy_tree)} #if not binary tree

# There is an issue with the names of the nodes in the tree, that is, some of the names are repeated even though they are from different branches in the tree. So, will need to change the names of the nodes to those that are unique

phy_tree(phylo) = ape::makeNodeLabel(phy_tree(phylo), method="number", prefix= 'n')

# Add 1 to avoid fractions on a zero denominator    

data.no0 = transform_sample_counts(phylo, function(x) x+1)
phylot = merge_phyloseq(data.no0,phylo@phy_tree,phylo@sam_data,phylo@tax_table) # make your new GP phyloseq object based on the newly created matrix

myMatrix = phylot@otu_table@.Data %>% t(.)
tree = phy_tree(phylot)

gp.philr = philr::philr(myMatrix, tree, part.weights='uniform', ilr.weights='uniform') %>% suppressMessages(.)
gp.dist = dist(gp.philr, method="euclidean") # Make it into phyloseq compatible object
}else{
  
  gp.dist = vegan::vegdist(t(mytsv),method=dissimilarity)  
  
}

#Write our the distribution between the different samples
write.table(as.matrix(gp.dist), sep = "\t", file = paste0(output),col.names=NA,quote=FALSE)

