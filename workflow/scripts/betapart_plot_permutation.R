
set.seed(1234)
suppressWarnings(suppressMessages(library(optparse)))	
suppressWarnings(suppressMessages(require(tidyverse)))	
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="folder with tsv file(s)", metavar="Folder with output table(s) from the permutations script"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Mapping file"),
  make_option(c("-c", "--color"), type="character", default=NULL, help="the suffix that is used for colors of the column", metavar="Color column suffix"),
  make_option(c("-b", "--basename"), type="character", default=NULL, help="Basename of the biom file", metavar="Basename of the biom file"),
  make_option(c("-d", "--distance"), type="character", default=NULL, help="Type of distance matrix. Either the string 'bray' or 'jaccard'.", metavar="Distance matrix"),
  make_option(c("-x", "--width"), type="character", default=NULL, help="Width for the graph", metavar="Graph width"),
  make_option(c("-y", "--height"), type="character", default=NULL, help="Height for the graph", metavar="Graph height"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder to save the files", metavar="Output folder")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
distance = opt$distance

myFolder = opt$input 
myFiles = list.files(path = myFolder)

output = opt$output
color = opt$color 

mapping = opt$mapping

mapping = suppressMessages(read.delim(mapping))
myGroups = c()

myheight = as.numeric(opt$height)
mywidth = as.numeric(opt$width)

basename = opt$basename


for(myfile in myFiles){
  s1 = unlist(strsplit(myfile, split=paste0("+"), fixed=TRUE))[2]
  myGroups = c(myGroups, unlist(strsplit(s1, split='–permutations.txt', fixed=TRUE))[1])
}

myGroups = unique(myGroups)

condition = unlist(strsplit(myFiles[1], split=paste0("+"), fixed=TRUE))[1]

concatLongTable = NULL

# Let's now filter the files based on the groups
for (thisGroup in myGroups) {

  mappingCategory = select(mapping, paste(condition))
  mappingColor    = select(mapping, paste0(condition,color))
  filteredMapping = cbind(mappingCategory,mappingColor) %>% unique(.)
  
  currentFiles = list.files(path = myFolder, pattern = thisGroup)
 
   currentFileTable = suppressMessages(readr::read_tsv(paste0(myFolder,condition,"+",thisGroup,"–permutations.txt"))) %>% as.data.frame(.)
   
   currentColor = which(filteredMapping[,1] == thisGroup) %>% filteredMapping[.,2] %>% as.character(.)
   currentCategory = rep(thisGroup, times = dim(currentFileTable)[1])
   
   if (distance == "bray"){
     colnames(currentFileTable) = c("Balanced","Gradient","BrayCurtis")
   }else{
     colnames(currentFileTable) = c("Turnover","Nestedness","Jaccard")
    }
   
   
   longTable = tidyr::pivot_longer(cbind(currentCategory,currentFileTable),
                                   cols= -currentCategory)
   currentColor = rep(currentColor, times = dim(longTable)[1])
   longTable = cbind(longTable,currentColor)
   
   
   colnames(longTable) = c("Category","Type","Value","Color")
     concatLongTable = rbind(concatLongTable,longTable)
 }
  

myPlot =  ggplot2::ggplot(concatLongTable,aes(x=Value,color=Category, linetype=Type))+
     geom_density()+
     theme_minimal()+
     ggtitle(paste0("Bray-Curtis Breakdown for ", thisGroup))+
     scale_color_manual(values = unique(as.character(concatLongTable$Color)))+
     scale_linetype_manual(values=c("twodash","solid", "dotted"))


ggsave(filename=paste0(output),plot=myPlot,width=mywidth,height=myheight,device="svg")
 

 
   