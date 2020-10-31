#!/usr/local/bin/Rscript --vanilla

set.seed(1234)

suppressWarnings(suppressMessages(library(phyloseq)))	
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(philr)))	
suppressWarnings(suppressMessages(require(ape)))
suppressWarnings(suppressMessages(library(optparse)))	

# Main function
main = function(phylot,myLevels,fileName){
  #phylot=phylo
  #myLevels=condName
  #fileName=output
  
  #====== The functions =======
  prunePhylo = function(phylot,myLevels){
    samNames = as.data.frame(sample_data(phylot))
    myNames = subset(samNames, myLevels %in% myLevels) %>% .$X.SampleID %>% as.character(.)
    mySubset = prune_samples(myNames,phylot)
    mySubset = prune_taxa(taxa_sums(mySubset) != 0, mySubset)  # Remove zero taxa
  }
  psuedoCount = function(phylot){

    data.no0 = transform_sample_counts(phylot, function(x) x+1)
    phylot = merge_phyloseq(data.no0,phylot@phy_tree,phylot@sam_data,phylot@tax_table) # make your new GP phyloseq object based on the newly created matrix
    
    return(phylot)
  }

    
  #===== The Steps ======
  
  # Prune your phyloseq object to only those with the requested input characteristics, and remove the taxa with zero values in all samples
  mySubset = prunePhylo(phylo,myLevels)

  # Add 1 to avoid fractions on a zero denominator
  myMatrix = psuedoCount(mySubset)
  myMatrix = myMatrix@otu_table@.Data %>% t(.)
  tree = phy_tree(phylot)
  
  gp.philr = philr(myMatrix, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
  gp.dist = dist(gp.philr, method="euclidean") # Make it into phyloseq compatible object
  
  #Write our the distribution between the different samples
  write.table(as.matrix(gp.dist), sep = "\t", file = paste0(fileName), na = "\t")

}

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Biom file", metavar="Features input file formatted as biom"),
  make_option(c("-t", "--tree"), type="character", default=NULL, help="Phylogenic Tree file. With or without lineage distances", metavar="Phylogenic tree file"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-g", "--groups"), type="character", default=NULL, help="Column name in the mapping file that indicates which grouping the samples belong to", metavar="Group indication column"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file name", metavar="Output file name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load files

biom = opt$input
biom = import_biom(biom)

map = opt$mapping 
map = import_qiime_sample_data(map)

tree = opt$tree 
tree = read_tree(tree)

condition = opt$groups
output = opt$output 

phylo = merge_phyloseq(biom,tree,map)

rm(biom,tree,map)

# Check if the tree is rooted and is binary (that is, divides only as a binary). If not, use the multi2di function from ape package to replace the multichotomies with a series of dichotomies with a branch length of zero.
if(!ape::is.rooted(phy_tree(phylo))){stop("Tree needs to be rooted.",call.=FALSE)} #if not rooted
if(!ape::is.binary.tree(phy_tree(phylo))){phy_tree(phylo) = multi2di(phylo@phy_tree)} #if not binary tree

# There is an issue with the names of the nodes in the tree, that is, some of the names are repeated even though they are from different branches in the tree. So, will need to change the names of the nodes to those that are unique
phy_tree(phylo) = ape::makeNodeLabel(phy_tree(phylo), method="number", prefix= 'n')

condNum = which(colnames(phylo@sam_data)==condition)
condName = phylo@sam_data[,condNum][[1]] %>% unique(.) %>% as.character(.)

nonZeroPhylo = suppressMessages(main(phylo,condName,output))	
