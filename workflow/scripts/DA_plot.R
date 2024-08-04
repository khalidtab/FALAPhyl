#!/usr/local/bin/Rscript --vanilla

library("optparse")
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Folder with the files to create the graph", metavar="Input folder"),
  make_option(c("-a", "--auc"), type="character", default=NULL, help="Path for AUC graph", metavar="Path for AUC graph"),
  make_option(c("-f", "--fdr"), type="character", default=NULL, help="Path for FDR graph", metavar="Path for FDR graph"),
  make_option(c("-p", "--power"), type="character", default=NULL, help="Path for Power graph", metavar="Path for Power graph"),
  make_option(c("-s", "--score"), type="character", default=NULL, help="Path for Score graph", metavar="Path for Scores graph"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output tsv file", metavar="Output tsv File name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

print(opt$input)

myFiles = list.dirs(path=opt$input) %>% 
          grep("effectSize",x=.,value = TRUE) %>% 
          sapply(., function(x) list.files(x, full.names = TRUE))

typesOfTests = gsub(pattern = ".*–[0-9]*", replacement = "",x=myFiles) %>% 
               gsub(pattern = "/.*", replacement = "",x=.) %>%
               gsub(pattern = "[0-9]–", replacement = "",x=.) %>%
               unique(.)
 

theTestTables = function(theTest){
  grep(theTest,myFiles,value = TRUE) %>%
               lapply(., function(x) read_tsv(x,show_col_types = FALSE)) %>%
               bind_rows(.) %>%
               as.data.frame(.)
}

theTestFiles = lapply(typesOfTests, function(x) theTestTables(x)) %>% bind_rows(.)

theTestFiles = theTestFiles %>%
  mutate(score = (AUC - 0.5) * Power - FDR)

write_tsv(theTestFiles,opt$output)

myAUC = ggerrorplot(theTestFiles, x = "Method", y = "AUC", 
            desc_stat = "mean", color = "Method" 
            , size=1,add = "jitter") + theme(axis.text.x = element_text(angle=90)) + labs(x="Method", y = "AUC") +
  geom_hline(yintercept=0.5,color="red") + rremove("legend")

myPower = ggerrorplot(theTestFiles, x = "Method", y = "Power", 
                    desc_stat = "mean", color = "Method" 
                    , size=1,add = "jitter") + theme(axis.text.x = element_text(angle=90)) + labs(x="Method", y = "Power")+ rremove("legend")

myFDR = ggerrorplot(theTestFiles, x = "Method", y = "FDR", 
                      desc_stat = "mean", color = "Method" 
                      , size=1,add = "jitter") + theme(axis.text.x = element_text(angle=90)) + labs(x="Method", y = "False Discovery Rate", caption="Higher AUC means spiked features have been identified,\nand AUC=0.5 means spiked features are randomly spread with non-spiked.\nTherefore, we want an AUC as high as possible.\nThat is, spiked features should have low p-values (ie they have been identified).\nPower is the proportion of spiked features that are significant after multiple p-value corrections\nIt is therefore the proportion of features you would expect to detect in a regular analysis.\nThe higher the power, the better.\nFDR indicates the proportion of significant features (after multiple correction)\nthat were not spiked and therefore shouldn't be significant.\nThis should be as low as possible.")+ rremove("legend")


myScores = ggerrorplot(theTestFiles, x = "Method", y = "score", 
                       desc_stat = "mean", color = "Method" 
                       , size=1,add = "jitter") + theme(axis.text.x = element_text(angle=90)) + labs(x="Method", y = "Score", title= 
                                                                                                       "Score", caption="Score is calculated for each method as follows: (Area Under the ROC Curve - 0.5) * Power - False Discovery Rate.\nThe higher the Score, the better the method is estimated to be.")+ rremove("legend")





svg(opt$auc,width = 10, height = 20)
myAUC
dev.off()

svg(opt$power,width = 10, height = 20)
myPower
dev.off()

svg(opt$fdr,width = 10, height = 20)
myFDR
dev.off()

svg(opt$score ,width = 10, height = 20)
myScores
dev.off()