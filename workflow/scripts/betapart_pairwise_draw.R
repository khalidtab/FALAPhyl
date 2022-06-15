#!/usr/local/bin/Rscript --vanilla

set.seed(1234)

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(require(tidyverse)))
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(effectsize)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Melted distance matrix", metavar="Bray-Curtis or Jaccard Distance melted matrix"),
  make_option(c("-d", "--distance"), type="character", default=NULL, help="Type of distance matrix. Either the string 'bray' or 'jaccard'.", metavar="Distance matrix"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="Mapping file", metavar="Input mapping file"),
  make_option(c("-c", "--category"), type="character", default=NULL, help="Category", metavar="Column in the mapping file that is the category to be compared"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output folder for all pairwise files", metavar="Output folder name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load variables
cond = opt$category
dist_type = opt$distance

# Load table
dist = opt$input %>% read_tsv(.,progress = FALSE)

currentConds = basename(opt$input) %>% str_split(.,".tsv")%>% .[[1]] %>% .[1] %>% str_split(.,"–") 
myCond1 = currentConds[[1]][1]
myCond2 = currentConds[[1]][2]

# Load mapping file
map = opt$mapping 
map = suppressMessages(read_tsv(map))

output = opt$output 

column1 = dist[,1] %>% as.data.frame(.) %>% .[,1] %>% as.numeric(.)
column2 = dist[,2] %>% as.data.frame(.) %>% .[,1] %>% as.numeric(.)

dist2 = cbind(rownames(dist),dist)
colnames(dist2)[1] = "SampleID"
dist3 = dist2 %>% pivot_longer(!SampleID,names_to = "type",values_to="value")

effectSize = rank_biserial(column1,column2, paired = TRUE) %>% interpret(., rules = "funder2019")
pvalue = wilcox.test(column1,column2, paired = TRUE) %>% .$p.value

if (pvalue == 0){ # R sets the pvalue to zero if it is below a certain level, so put it as at least 0.0001 if it is too low
  pvalue = 0.0001
}

if (dist_type == "bray"){
  myGraph = ggpubr::ggpaired(dist, cond1 = "balanced variation in abundance", cond2 = "abundance gradients", fill = "condition", palette = "jco", line.size=0.01, 
                             title= (paste0(dist_type, " breakdown\n",myCond1," and ",myCond2, "\np-value =",formatC(pvalue, format = "g", digits = 2),  ", Rank biserial effect size=",formatC(effectSize$r_rank_biserial, format = "g", digits = 2))))
} else {
  myGraph = ggpubr::ggpaired(dist, cond1 = "Turn-over", cond2 = "Nestedness", fill = "condition", palette = "jco", line.size=0.01, 
                             title= (paste0(dist_type, " breakdown\n",myCond1," and ",myCond2, "\np-value =",formatC(pvalue, format = "g", digits = 2),". Rank biserial effect size= ",formatC(effectSize$r_rank_biserial, format = "g", digits = 2))))  
} 


ggsave(filename=paste0(output),plot=myGraph)


