# matrix to pairwise
library("optparse")
#library(reshape2)
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="dissimilarity matrix", metavar="Input dissimilarity matrix"),
  make_option(c("-p", "--pvalue"), type="character", default=NULL, help="p-value matrix", metavar="pvalue file"),
  make_option(c("-n", "--nodes"), type="character", default=NULL, help="Name of the output nodes table", metavar="output nodes table file name"),
  make_option(c("-e", "--edges"), type="character", default=NULL, help="Name of the output edges table", metavar="output edges table file name"),
  make_option(c("-a", "--threshold"), type="character", default=NULL, help="Correlation cut off threshold below which correlations will be discarded", metavar="Correlation threshold"),  
  make_option(c("-b", "--pvaluethreshold"), type="character", default=NULL, help="P-value cut off threshold below which pvalues will be discarded", metavar="P-value threshold")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

input = opt$input #"/Users/khaled/Desktop/bioinfo_snakemake/data/network/kfp5s2_CPd_GAPd_LAPd_corr/corr+CPd.tsv"
pvalue = opt$pvalue #"/Users/khaled/Desktop/bioinfo_snakemake/data/network/kfp5s2_CPd_GAPd_LAPd_corr/pvalue+CPd.tsv"
threshold = opt$threshold #"0.6"
threshold = as.numeric(threshold)
pvaluethreshold = opt$pvaluethreshold #"0.1"
pvaluethreshold = as.numeric(pvaluethreshold)

nodesOutput = opt$nodes #"/Users/khaled/Desktop/bioinfo_snakemake/data/network/kfp5s2_CPd_GAPd_LAPd_corr/nodes+CPd.tsv"
edgesOutput = opt$edges #"/Users/khaled/Desktop/bioinfo_snakemake/data/network/kfp5s2_CPd_GAPd_LAPd_corr/edges+CPd.tsv"

# Unpack/Melt the correlation matrix table
myInput = suppressMessages(read_tsv(input))
myInput =  as.data.frame(myInput)
row.names(myInput) = myInput[,1]
myInput$`#OTU ID` = NULL # by now, row names and column names are the feature ID

myInput = as.matrix(myInput)
myInput = replace(myInput, lower.tri(myInput, TRUE), NA)
myInput = reshape2::melt(myInput, na.rm = TRUE)
colnames(myInput) = c("Source","Target","corrWeight")

# Unpack/Melt the pvalue matrix table
mypvalue = suppressMessages(read_tsv(pvalue))
mypvalue =  as.data.frame(mypvalue)
row.names(mypvalue) = mypvalue[,1]
mypvalue$`#OTU ID` = NULL # by now, row names and column names are the feature ID

mypvalue = as.matrix(mypvalue)
mypvalue = replace(mypvalue, lower.tri(mypvalue, TRUE), NA)
mypvalue = reshape2::melt(mypvalue, na.rm = TRUE)
colnames(mypvalue) = c("Source","Target","pvalue")

myTable = merge(myInput,mypvalue)
Type = rep("Undirected", length(rownames(myTable)))
myTable = cbind(myTable,Type)

myTable = cbind(myTable, abs(myTable$corrWeight))
myTable = subset(myTable, abs(myTable$corrWeight) >= threshold)
myTable$`abs(myTable$corrWeight)` = NULL

weight = myTable$corrWeight
myTable = cbind(myTable,weight)

myTable$weight[myTable$weight < 0] = 0 # This way, the louvain modularity function below will not consider negative correlations between two edges as a positive correlation and be influenced to place them in the same cluster

row.names(myTable) = NULL # At this point, you have a pairwise table that was filtered to the level specified

myigraph = graph.data.frame(myTable, directed = FALSE)
myModularity = cluster_louvain(myigraph) #Same method as Gephi

nodesList = cbind(myModularity$names,myModularity$membership)
colnames(nodesList) = c("Id","Community")
nodesList = as.data.frame(nodesList)

# To make it easier to do the Zi-Pi calculations for each OTU, copy myTable to another, and append it with the "Target" OTU being the "Source"
myTable$weight = NULL
newTable = myTable
colnames(newTable) = c("Target","Source","corrWeight","pvalue","Type")
newTable = cbind(newTable$Source,newTable$Target,
                 newTable$corrWeight,newTable$pvalue,
                 as.character(newTable$Type))
colnames(newTable) = c("Source","Target","corrWeight","pvalue","Type")
newTable = as.data.frame(newTable, stringsAsFactors=FALSE)

myTable2 = rbind(newTable,myTable)

#Combine nodesList to the myTable2 to make it easier subset the table
Source = nodesList %>% as.data.frame()
colnames(Source) = c("Source","SourceCommunity")
myTable2 = merge(myTable2,Source)
Target = nodesList %>% as.data.frame()
colnames(Target) = c("Target","TargetCommunity")
myTable2 = merge(myTable2,Target)
colnames(myTable2) = c("Target","Source","Weight","pvalue",
                       "Type","SourceCommunity","TargetCommunity")
myTable = myTable2
# Before Zi-Pi plot calculations, a little variable clean up
rm(myTable2,Source, Target, newTable,myigraph,
   threshold, pvaluethreshold,Type,mypvalue,myInput,weight,input,pvalue)

# Create Zi-Pi table
final_table = matrix(NA, ncol = 3, nrow = length(nodesList$Id)) %>% as.data.frame(.,stringsAsFactors=FALSE)
colnames(final_table) = c("Id","Zi","Pi")
final_table$Id = nodesList$Id


myCommunities = unique(myTable$SourceCommunity) %>% as.character(.)

#Zi calculations
for (x in (1:length(myCommunities))){

## Leave only edges coming nodes from the same community
myComm = myCommunities[x] %>% as.character(.)
currentComm = subset(myTable, SourceCommunity==myComm)
currentComm = subset(currentComm, TargetCommunity==myComm)

if (dim(currentComm)[1] == 0){next} #Break out if the current Community has only negative edges to other nodes, or if the node doesn't have any edges to nodes within the same community

uniq_comm_ids = unique(currentComm$Source)
comm_zi = matrix(NA, ncol = 5, nrow = length(uniq_comm_ids)) %>% as.data.frame(.,stringsAsFactors=FALSE)
colnames(comm_zi) = c("Id","NumOfConnections","average","STDev","Zi")
comm_zi$Id = uniq_comm_ids

for (xx in uniq_comm_ids){

  myNum = subset(currentComm, Source==xx) 
  myNum = dim(myNum)
  myNum = myNum[1]
  comm_zi = within(comm_zi, NumOfConnections[Id == xx] <- myNum)
}

comm_zi_STDEV = sd(comm_zi$NumOfConnections)

if (comm_zi_STDEV == 0){
  comm_zi$Zi = NA
} else {
  comm_zi$STDev = comm_zi_STDEV
  comm_zi$average = ave(comm_zi$NumOfConnections)
  comm_zi$Zi = (comm_zi$NumOfConnections - comm_zi$average)/comm_zi$STDev
}

comm_zi$NumOfConnections = NULL
comm_zi$average = NULL
comm_zi$STDev = NULL

final_table = merge(final_table,comm_zi, by="Id", all.x = TRUE, all.y = TRUE)
myValues = coalesce(as.double(final_table$Zi.y),as.double(final_table$Zi.x))
final_table$Zi.x = NULL
colnames(final_table) = c("Id","Pi","Zi")
final_table$Zi = myValues
}

#Clean up after Zi calculations
rm(comm_zi,comm_zi_STDEV,currentComm,myComm,myCommunities,myNum,myValues,uniq_comm_ids,x,xx)

#Pi calculations
for (x in nodesList$Id){
  piTotalTable = subset(myTable, Source==x)
  total = dim(piTotalTable)[1]
  sigma = 0
  
  for (y in unique(piTotalTable$TargetCommunity)){
    commConnections = subset(piTotalTable, TargetCommunity==y) %>% dim(.) %>% .[1]
    fraction = (commConnections/total)^2
    sigma = sigma + fraction
  }
  
  featurePiScore = 1-sigma
  final_table = within(final_table, Pi[Id == x] <- featurePiScore)
} 

final_table = merge(final_table,nodesList)

# final_table is nodes table, myTable is edges table

rm(x,y,total,sigma,piTotalTable,myModularity,fraction,featurePiScore,nodesList,commConnections)

write.table(final_table, file = nodesOutput,sep = "\t", row.names = FALSE)
write.table(myTable, file = edgesOutput,sep = "\t", row.names = FALSE)

