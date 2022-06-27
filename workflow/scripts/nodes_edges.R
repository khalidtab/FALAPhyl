# matrix to pairwise
library("optparse")
#library(reshape2)
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(igraph)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="correlation matrix", metavar="Input correlation matrix"),
  make_option(c("-p", "--pvalue"), type="character", default=NULL, help="p-value matrix", metavar="pvalue file"),
  make_option(c("-n", "--nodes"), type="character", default=NULL, help="Name of the output nodes table", metavar="output nodes table file name"),
  make_option(c("-e", "--edges"), type="character", default=NULL, help="Name of the output edges table", metavar="output edges table file name"),
  make_option(c("-a", "--threshold"), type="character", default=NULL, help="Correlation cut off threshold below which correlations will be discarded", metavar="Correlation threshold"),  
  make_option(c("-b", "--pvaluethreshold"), type="character", default=NULL, help="P-value cut off threshold above which pvalues will be discarded", metavar="P-value threshold")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

input = opt$input
pvalue = opt$pvalue
threshold = opt$threshold
threshold = as.numeric(threshold)
pvaluethreshold = opt$pvaluethreshold
pvaluethreshold = as.numeric(pvaluethreshold)

nodesOutput = opt$nodes
edgesOutput = opt$edges

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
myTable = subset(myTable, abs(myTable$corrWeight) >= threshold) # Subset to only correlations that fit the criteria
myTable$`abs(myTable$corrWeight)` = NULL

myTable = subset(myTable, abs(myTable$pvalue) < pvaluethreshold) # Subset to only pvalues less than the pvalue threshold

weight = myTable$corrWeight
myTable = cbind(myTable,weight)

myTable$weight[myTable$weight < 0] = 0 # This way, the louvain modularity function below will not consider negative correlations between two edges as a positive correlation and be influenced to place them in the same cluster
row.names(myTable) = NULL # At this point, you have a pairwise table that was filtered to the level specified

tableForModularity = myTable
tableForModularity = subset(tableForModularity, tableForModularity$weight != 0) # Subset to only correlations that are positive
tableForModularity = subset(tableForModularity, tableForModularity$Source != "Others") # Remove others so it doesn't mess up the Louvain algorithm

myigraph = graph.data.frame(tableForModularity, directed = FALSE)
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
   Type,mypvalue,myInput,weight,input,pvalue)

# Create Zi-Pi table
zipiTable = matrix(NA, ncol = 3, nrow = length(nodesList$Id)) %>% as.data.frame(.,stringsAsFactors=FALSE)
colnames(zipiTable) = c("Id","Zi","Pi")
zipiTable$Id = nodesList$Id


myCommunities = unique(myTable$SourceCommunity) %>% as.character(.)

#Zi calculations
for (x in (1:length(myCommunities))){

## Leave only edges coming from nodes in the same community
myComm = myCommunities[x] %>% as.character(.)
currentComm = subset(myTable, SourceCommunity==myComm)
currentComm = subset(currentComm, TargetCommunity==myComm)

if (dim(currentComm)[1] == 0){next} #Break out if the current Community has only negative edges to other nodes, or if the node doesn't have any edges to nodes within the same community. These communities don't get reported through this method.


uniq_comm_ids = c(unique(currentComm$Source),unique(currentComm$Target))
comm_zi = matrix(NA, ncol = 7, nrow = length(uniq_comm_ids)) %>% as.data.frame(.,stringsAsFactors=FALSE)
colnames(comm_zi) = c("Id","NumOfSourceConnections","NumOfTargetConnections","NumOfConnections","average","STDev","Zi")
comm_zi$Id = uniq_comm_ids

for (xx in uniq_comm_ids){ # For source

  myNum = subset(currentComm, Source==xx) 
  myNum = dim(myNum)
  myNum = myNum[1]
  comm_zi = within(comm_zi, NumOfSourceConnections[Id == xx] <- myNum)
}

for (xx in uniq_comm_ids){ # For target
  
  myNum = subset(currentComm, Target==xx) 
  myNum = dim(myNum)
  myNum = myNum[1]
  comm_zi = within(comm_zi, NumOfTargetConnections[Id == xx] <- myNum)
}

comm_zi$NumOfConnections = comm_zi$NumOfSourceConnections + comm_zi$NumOfTargetConnections
comm_zi$NumOfSourceConnections = NULL
comm_zi$NumOfTargetConnections = NULL

comm_zi_STDEV = sd(comm_zi$NumOfConnections)


if (is.na(comm_zi_STDEV) == TRUE){
  comm_zi$Zi = 0 # Divide by zero is not allowed, therefore we are assigning it as zero instead.
  } else if(comm_zi_STDEV == 0) {
    comm_zi$Zi = 0 # Divide by zero is not allowed, therefore we are assigning it as zero instead.
  } else {
  comm_zi$STDev = comm_zi_STDEV
  comm_zi$average = ave(comm_zi$NumOfConnections)
  comm_zi$Zi = (comm_zi$NumOfConnections - comm_zi$average)/comm_zi$STDev
}

comm_zi$NumOfConnections = NULL
comm_zi$average = NULL
comm_zi$STDev = NULL

zipiTable = merge(zipiTable,comm_zi, by="Id", all.x = TRUE, all.y = TRUE)
myValues = coalesce(as.double(zipiTable$Zi.y),as.double(zipiTable$Zi.x))
zipiTable$Zi.x = NULL
colnames(zipiTable) = c("Id","Pi","Zi")
zipiTable$Zi = myValues
}

#Clean up after Zi calculations
rm(comm_zi,comm_zi_STDEV,currentComm,myComm,myCommunities,myNum,myValues,uniq_comm_ids,x,xx)

zipiTable = unique(zipiTable) # There are some duplicate measurements of the table, so only keep unique ones

#Pi calculations
for (x in nodesList$Id){
  # x = nodesList$Id[2]
  piTotalTableSource = subset(myTable, Source==x)
  piTotalTableTarget = subset(myTable, Target==x)
  colnames(piTotalTableTarget) = c("Source","Target","Weight","pvalue","Type","TargetCommunity","SourceCommunity") # Change the target to a source, including the community
  
  piTotalTable = rbind(piTotalTableSource,piTotalTableTarget)
  
  total = dim(piTotalTable)[1]
  sigma = 0
  
  for (y in unique(piTotalTable$TargetCommunity)){
    commConnections = subset(piTotalTable, TargetCommunity==y) %>% dim(.) %>% .[1]
    fraction = (commConnections/total)^2
    sigma = sigma + fraction
  }
  
  featurePiScore = 1-sigma
  zipiTable = within(zipiTable, Pi[Id == x] <- featurePiScore)
} 

zipiTable = merge(zipiTable,nodesList)
# zipiTable2 = dplyr::full_join(zipiTable,nodesList)
# zipiTable is nodes table, myTable is edges table

rm(x,y,total,sigma,piTotalTable,myModularity,fraction,featurePiScore,nodesList,commConnections)

write.table(zipiTable, file = nodesOutput,sep = "\t", row.names = FALSE)
write.table(myTable, file = edgesOutput,sep = "\t", row.names = FALSE)


