
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
myFiles = list.files(path = myFolder, pattern = paste0("^betapart_permutations_",distance,"_"))

output = opt$output
color = opt$color 

mapping = opt$mapping

mapping = suppressMessages(read.delim(mapping))
myGroups = c()

basename = opt$basename


for(myfile in myFiles){
  s1 = unlist(strsplit(myfile, split=paste0('betapart_permutations_',distance,"_"), fixed=TRUE))[2]
  myGroups = c(myGroups, unlist(strsplit(s1, split='+', fixed=TRUE))[1])  
}

myGroups = unique(myGroups)

#Functions
getCurrentFileTable = function(currentFile){
  
  currentCategory = unlist(strsplit(currentFile, split='betapart_permutations_', fixed=TRUE))[2]
  currentCategory = unlist(strsplit(currentCategory, split='+', fixed=TRUE))[2]
  currentCategory = unlist(strsplit(currentCategory, split='.tsv', fixed=TRUE))[1]
  
  currentFileTable = suppressMessages(readr::read_tsv(paste0(myFolder,currentFile))) %>% as.data.frame(.)
  colnames(currentFileTable) = c("Balanced","Gradient","BrayCurtis")
  
  currentColor = which(filteredMapping[,1] == currentCategory) %>% filteredMapping[.,2] %>% as.character(.)
  currentCategory = rep(currentCategory, times = dim(currentFileTable)[1])
  
  longTable = tidyr::pivot_longer(cbind(currentCategory,currentFileTable),
                                  cols= -currentCategory)
  currentColor = rep(currentColor, times = dim(longTable)[1])
  longTableColor = cbind(longTable,currentColor)
  
  
  
  
  return(longTableColor)
}
#########


# Let's now filter the files based on the groups
for (thisGroup in myGroups) {

  mappingCategory = select(mapping, paste(thisGroup))
  mappingColor    = select(mapping, paste0(thisGroup,color))
  filteredMapping = cbind(mappingCategory,mappingColor) %>% unique(.)
  
  currentFiles = list.files(path = myFolder, pattern = thisGroup)
  concatLongTable = NULL
 
 for (currentFile in currentFiles) {
   currentLongTable = getCurrentFileTable(currentFile)
   colnames(currentLongTable) = c("Category","Type","Value","Color")
     concatLongTable = rbind(concatLongTable,currentLongTable)
 }
  

myPlot =  ggplot2::ggplot(concatLongTable,aes(x=Value,color=Category, linetype=Type))+
     geom_density()+
     theme_minimal()+
     ggtitle(paste0("Bray-Curtis Breakdown for ", thisGroup))+
     scale_color_manual(values = unique(as.character(concatLongTable$Color)))+
     scale_linetype_manual(values=c("twodash","solid", "dotted"))


ggsave(filename=paste0(output,"betapart_",basename,"+",thisGroup,".svg"),plot=myPlot)
   
   }
 

 
   