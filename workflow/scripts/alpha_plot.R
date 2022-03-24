
library("optparse")
suppressMessages(library("vegan"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("cowplot"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="alpha div output", metavar="Input alpha div output"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output SVG file", metavar="Output SVG File name"),
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

alpha = opt$input
output = opt$output
category = opt$group
mapping_file = opt$mapping 
color = opt$color

alpha_table = suppressMessages(read_tsv(alpha))
typeOfAlpha = colnames(alpha_table)[2]
map = suppressMessages(read.csv(mapping_file, skip=0, header=T, sep="\t"))


test = merge(alpha_table,map,by.x="Samples",by.y="X.SampleID",sort = TRUE)

colNames = test[which( colnames(test)==category)] %>% as.matrix(.) %>% as.character(.) %>% unique(.) 

colNum = which( colnames(test)==category) #for category
condNum = which( colnames(test)==typeOfAlpha) #for alpha measurements
colorNum = which( colnames(test)==color) #for color column

if (length(colNames) == 2){
mycomparisons = as.matrix(colNames) %>% as.character(.)
mycomparisons = list(mycomparisons)
} else
{
mycomparisons = combn(colNames,m=2,simplify=FALSE)  
}

mycomparisons1st = NULL
mycomparisons2nd = NULL

for (i in mycomparisons){
group1=subset(test,test[colNum] == i[1]) %>% .[condNum] %>% as.matrix(.) %>% as.numeric(.)
group2=subset(test,test[colNum] == i[2]) %>% .[condNum] %>% as.matrix(.) %>% as.numeric(.)

wilcoxtest=wilcox.test(group1,group2,exact=FALSE)

if (wilcoxtest$p.value < 0.05){
  mycomparisons1st = list(mycomparisons1st, i[1])
  mycomparisons2nd = list(mycomparisons2nd, i[2])
}}

mycomparisons1st = unlist(mycomparisons1st)
mycomparisons2nd = unlist(mycomparisons2nd)


mycomparisonsCurated = list()
if (length(mycomparisons1st) > 0){
for (i in 1:length(mycomparisons1st)){ 
  mycomparisonsCurated[[i]] = c(mycomparisons1st[i],mycomparisons2nd[i])
}
}

mycolors = test[colorNum] %>% as.matrix(.) %>% as.character(.) %>% unique(.)

myplot = ggviolin(test, x=category,y=typeOfAlpha,fill=category,rug=TRUE,color=category,
          alpha=0.5,palette = mycolors,
         add = "boxplot", add.params = list(fill = "white"), 
         )+
    stat_compare_means(comparisons=mycomparisonsCurated,
                       method="wilcox",hide.ns = TRUE) +
  theme(legend.position = "none") + xlab(NULL)

ggsave(output,plot=myplot)
