#!/usr/local/bin/Rscript --vanilla

library("optparse")
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Folder with the files to create the graph", metavar="Input folder"),
  make_option(c("-g", "--graph"), type="character", default=NULL, help="Name of the output SVG file", metavar="Output SVG File name"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output tsv file", metavar="Output tsv File name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

myFiles = list.files(path=opt$input,pattern="diffPower–")

readFiles = data.frame(matrix(ncol = 6, nrow = 0))

for (myFile in myFiles){
  myReadFile = suppressMessages(read_tsv(paste0(opt$input,myFile)))
  readFiles = rbind(hi,readFiles)
}

myReadFile = readFiles %>%
  mutate(score = (AUC - 0.5) * Power - FDR)

write_tsv(myReadFile,opt$output)

myAUC = ggerrorplot(myReadFile, x = "Method", y = "AUC", 
            desc_stat = "mean", color = "Method" 
            , size=1,add = "jitter") + theme(axis.text.x = element_text(angle=90)) + labs(x="Method", y = "AUC") +
  geom_hline(yintercept=0.5,color="red") + rremove("legend")

myPower = ggerrorplot(myReadFile, x = "Method", y = "Power", 
                    desc_stat = "mean", color = "Method" 
                    , size=1,add = "jitter") + theme(axis.text.x = element_text(angle=90)) + labs(x="Method", y = "Power")+ rremove("legend")

myFDR = ggerrorplot(myReadFile, x = "Method", y = "FDR", 
                      desc_stat = "mean", color = "Method" 
                      , size=1,add = "jitter") + theme(axis.text.x = element_text(angle=90)) + labs(x="Method", y = "False Discovery Rate", caption="Higher AUC means spiked features have been identified,\nand AUC=0.5 means spiked features are randomly spread with non-spiked.\nTherefore, we want an AUC as high as possible.\nThat is, spiked features should have low p-values (ie they have been identified).\nPower is the proportion of spiked features that are significant after multiple p-value corrections\nIt is therefore the proportion of features you would expect to detect in a regular analysis.\nThe higher the power, the better.\nFDR indicates the proportion of significant features (after multiple correction)\nthat were not spiked and therefore shouldn't be significant.\nThis should be as low as possible.")+ rremove("legend")


myScores = ggerrorplot(myReadFile, x = "Method", y = "score", 
                       desc_stat = "mean", color = "Method" 
                       , size=1,add = "jitter") + theme(axis.text.x = element_text(angle=90)) + labs(x="Method", y = "Score", title= 
                                                                                                       "Score", caption="Score is calculated for each method as follows: (Area Under the ROC Curve - 0.5) * Power - False Discovery Rate.\nThe higher the Score, the better the method is estimated to be.")+ rremove("legend")



figure = ggarrange(myAUC, myPower, myFDR,myScores, ncol = 1, nrow = 4)

svg(opt$graph,width = 10, height = 20)
figure
dev.off()

