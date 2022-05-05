library("optparse")
suppressMessages(library("vegan"))
suppressMessages(library("tidyverse"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="dissimilarity matrix", metavar="Input dissimilarity matrix"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output for test", metavar="Output for test"),
  make_option(c("-b", "--boxplot"), type="character", default=NULL, help="Output for boxplot", metavar="Output for boxplot"),
  make_option(c("-p", "--pcoa"), type="character", default=NULL, help="Output for PCoA of beta dispersion with standard deviation ellipses", metavar="Beta dispersion PCoA"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="The mapping file", metavar="Mapping file"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="The category in the mapping file", metavar="Group name"),
  make_option(c("-c", "--color"), type="character", default=NULL, help="color suffix", metavar="color suffix"),
  make_option(c("-x", "--width"), type="character", default=NULL, help="Width of the SVG", metavar="Width of SVG"),
  make_option(c("-y", "--height"), type="character", default=NULL, help="Height of the SVG", metavar="Height of SVG"),
  make_option(c("-t", "--test"), type="character", default=NULL, help="Test to be done, options: adonis, betadisper", metavar="Output of test")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

dissimilarity = opt$input
mytest = opt$test
output = opt$output
category = opt$group
mapping_file = opt$mapping 
color = opt$color
output_for_boxplot = opt$boxplot
output_for_pcoa = opt$pcoa
mywidth = as.numeric(opt$width)
myheight = as.numeric(opt$height)

dis = read.table(dissimilarity,sep="\t",head=TRUE,row.names = 1) %>% as.dist(.)

# Prepare the mapping file
map = suppressMessages(read.csv(mapping_file, skip=0, header=T, sep="\t"))
colnames(map)[1] = "SampleID"
map = map[match(colnames(dis %>% as.matrix(.)), map$SampleID),] # Arrange your mapping file to have the same arrangement as your dissimilarity matrix


#Find what category the samples come from
catNum = which(colnames(map) == category)
colNum = which(colnames(map) == paste0(category,color))
working_map = cbind(as.character(map[,1]),
                    as.character(map[,catNum]),
                    as.character(map[,colNum])) %>% as.data.frame(.)
colnames(working_map) = c("SampleID","condition","color")



if (mytest == "betadisper"){
print("Beta-dispersion calculation starting…")
betadisper_test = vegan::betadisper(dis, working_map$condition, bias.adjust = TRUE) #Bias.adjust for small sample sizes is TRUE
print("Beta-dispersion permutation testing starting…")
betadisper_perm = betadisper_test %>% permutest(.,permutations=999,pairwise=TRUE)
print("Beta-dispersion Tukey HSD starting…")
pergroup        = betadisper_test %>% TukeyHSD(.,which="group",conf.level = 0.95)

arrangedColors = unique(working_map[ , c("condition", "color")]) %>% arrange(condition) %>% .$color # Get the correct order of colors based on the alphabetical order of the condition
arrangedCateg = unique(working_map[ , c("condition", "color")]) %>% arrange(condition) %>% .$condition # Get the correct order of categories based on the alphabetical order of the condition



svg(filename = paste0(output_for_boxplot),family="sans")
myboxplot = boxplot(betadisper_test,
                    notch=TRUE,outline=TRUE,
                    col=as.character(arrangedColors), 
                    xlab = paste0(category),
                    main = "Beta-dispersion distance from centroid",show.names=FALSE)
legend("topright", legend = myboxplot$names, border="black", fill = as.character(arrangedColors),cex=0.5)
dev.off()



svg(filename = paste0(output_for_pcoa),width = mywidth,height=myheight,family="sans")
p1 = plot(betadisper_test,hull=FALSE,label=FALSE,col=as.character(arrangedColors),
     ellipse=TRUE,lty="solid",lwd=5, 
     segments=TRUE,seg.col="grey",seg.lty=1,seg.lwd=0.3,
     pch=rep(16,length(arrangedColors)), # type of point
     xlab="PCoA1",ylab="PCoA2",
     main=paste("PCoA with spatial median and beta-dispersion standard deviation ellipses"))
legend("topright", legend = myboxplot$names, border="black", fill = as.character(arrangedColors),cex=0.5)
dev.off()

output_betadisper_perm = paste0(output)

suppressMessages(write.table(betadisper_perm$tab,output_betadisper_perm,
            sep="\t",append = TRUE,
            row.names = TRUE, 
            col.names = TRUE))

suppressMessages(write.table((betadisper_perm$pairwise %>% as.data.frame(.)),
            output_betadisper_perm,
            sep="\t",append = TRUE,
            row.names = TRUE, 
            col.names = TRUE))

} else if (mytest == "adonis"){
output_adonis = paste0(output)
output_pairwise_adonis = paste0(output,"_pairwise.txt")
print("ADONIS starting…")
myAdonis = vegan::adonis2(dis~condition, data=working_map)

myCombinations=combn(unique(working_map$condition),2)
pairwise_p = numeric()
pairwise_R2 = numeric()


for (i in 1:dim(myCombinations)[2]){
test1 = myCombinations[,i][1]
test2 = myCombinations[,i][2]


subset_map = working_map %>% filter(condition == test1 | condition == test2) 
combination_samples = subset_map %>% .[1] %>% unlist(.)
YesorNo = rownames(as.matrix(dis)) %in% combination_samples
subset_data = as.matrix(dis)[YesorNo,YesorNo] %>% as.dist(.)
subset_map = subset_map[match(colnames(subset_data %>% as.matrix(.)), subset_map$SampleID),] # Arrange your mapping file to have the same arrangement as your dissimilarity matrix

print(paste("ADONIS pairwise testing starting between",test1,test2))
subset_adonis = vegan::adonis2(subset_data~condition, data=subset_map)
pairwise_p[paste0(test1,"_",test2)]  = subset_adonis$`Pr(>F)`[1]
pairwise_R2[paste0(test1,"_",test2)]  = subset_adonis$R2[1]
}

pairwise_fdr_padj = p.adjust(pairwise_p,method="fdr")
pairwise_adonis_results = cbind(pairwise_R2,pairwise_p,pairwise_fdr_padj) %>% as.data.frame(.)

write.table(myAdonis,output_adonis,sep="\t",append = FALSE,row.names = TRUE, col.names = TRUE,quote=FALSE)
write.table(pairwise_adonis_results,output_pairwise_adonis,sep="\t",append = FALSE,row.names = TRUE, col.names = TRUE, quote=FALSE)

} else if (mytest == "anosim"){
  print("ANOSIM starting…")
  myANOSIM = vegan::anosim(dis, working_map$condition, permutations = 999)
  
  myCombinations=combn(unique(working_map$condition),2)
  pairwise_p = numeric()
  pairwise_R = numeric()
  
  for (i in 1:dim(myCombinations)[2]){
    test1 = myCombinations[,i][1]
    test2 = myCombinations[,i][2]
    
    subset_map = working_map %>% filter(condition == test1 | condition == test2) 
    combination_samples = subset_map %>% .[1] %>% unlist(.)
    YesorNo = rownames(as.matrix(dis)) %in% combination_samples
    subset_data = as.matrix(dis)[YesorNo,YesorNo] %>% as.dist(.)
    subset_map = subset_map[match(colnames(subset_data %>% as.matrix(.)), subset_map$SampleID),] # Arrange your mapping file to have the same arrangement as your dissimilarity matrix
    
    print(paste("ANOSIM pairwise combinations starting between",test1,test2))
    subset_anosim = vegan::anosim(subset_data, subset_map$condition, permutations = 999)
    pairwise_p[paste0(test1,"_",test2)]  = subset_anosim$signif
    pairwise_R[paste0(test1,"_",test2)]  = subset_anosim$statistic
  }
  pairwise_fdr_padj = p.adjust(pairwise_p,method="fdr")
  pairwise_anosim_results = cbind(pairwise_R,pairwise_p,pairwise_fdr_padj) %>% as.data.frame(.)
  
  output_anosim = paste0(output)
  output_pairwise_anosim = paste0(tools::file_path_sans_ext(output),"_pairwise.txt")
  
  sink(output_anosim)
  print(myANOSIM)
  sink()
  
  write.table(pairwise_anosim_results,output_pairwise_anosim,sep="\t",append = FALSE,row.names = TRUE, col.names = TRUE, quote=FALSE)
  
}
   