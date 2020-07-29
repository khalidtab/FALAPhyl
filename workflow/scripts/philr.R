#!/usr/bin/env Rscript --vanilla

library(optparse)
library(philr)
library(phyloseq)
library(ape)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="biom file", metavar="Input biom file"),
  make_option(c("-t", "--tree"), type="character", default=NULL, help="The tree file", metavar="Input tree file"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output dissimilarity matrix", metavar="Output dissimilarity matrix"),
  make_option(c("-v", "--viz"), type="character", default=NULL, help="Name of the output PCA in SVG format", metavar="Output SVG file"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="The mapping file", metavar="Mapping file"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="The category in the mapping file", metavar="Group name"),
  make_option(c("-c", "--color"), type="character", default=NULL, help="The color column in the mapping file", metavar="Color name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}


biom = opt$input # biom = "/Users/khaled/Desktop/bioinfo_snakemake/data/biom/Func_GAPd_LAPd_CPd_HNS_subsystem.biom"
tree = opt$tree # tree = "/Users/khaled/Desktop/bioinfo_snakemake/data/tree/Func_GAPd_LAPd_CPd_HNS_subsystem.tre"
output = opt$output # output = "/Users/khaled/Desktop/bioinfo_snakemake/data/distances/Philr_Func_GAPd_LAPd_CPd_HNS_subsystem.txt"
viz = opt$viz # viz = "/Users/khaled/Desktop/bioinfo_snakemake/data/plots/Philr_Func_GAPd_LAPd_CPd_HNS_subsystem.svg"
map = opt$mapping # map = "/Users/khaled/Desktop/bioinfo_snakemake/data/map/Func_GAPd_LAPd_CPd_HNS_subsystem.txt"
category = opt$group #category = "ShortCond"
color = opt$color #color = "Colors"


biom = import_biom(biom)
tree = read_tree(tree)
map = import_qiime_sample_data(map)
phylo = merge_phyloseq(biom,tree,map)
  
# make your new phyloseq object with taxa as columns
data = phylo@otu_table@.Data
data = t(data)
data = otu_table(data, taxa_are_rows = FALSE)
phylo = merge_phyloseq(data,phylo@phy_tree,phylo@sam_data,phylo@tax_table) 
  
# Add pseudocount of 0.5 to all counts so we're not dividing by zero
phylo = transform_sample_counts(phylo, function(x) x+0.5)
otu.table = (otu_table(phylo))
tree = phy_tree(phylo)
  
## Transform Data using PhILR
##Calculate the distances and (optionally) weight the resulting PhILR space using phylogenetic distance.
gp.philr = philr(otu.table, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
gp.dist = dist(gp.philr, method="euclidean") # Make it into phyloseq compatible object
  
#Write our the distribution between the different samples
test = as.data.frame(as.matrix(gp.dist))
test = cbind(rownames(test),test)
colnames(test)[1] = ""
write_tsv(test, path = output)


