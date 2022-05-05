library("optparse")
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(ggrepel)))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Nodes list with Zi-Pi coordinates", metavar="Zi-Pi input"),
  make_option(c("-t", "--threshold"), type="character", default=NULL, help="Correlation thershold", metavar="Correlation threshold"),
  make_option(c("-p", "--pvalue"), type="character", default=NULL, help="P-value", metavar="P-value"),
  make_option(c("-c", "--core"), type="character", default=NULL, help="The percentage (presented as a fraction of one, eg: 80% is 0.8) by which the feature needs to be present in a group to be considered as part of the core.", metavar="Core"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Path to output svg graph", metavar="Path of Output SVG")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

zipiTable = opt$input 
zipiTable = suppressMessages(read.csv(zipiTable, skip=0, header=T, sep="\t"))

threshold = opt$threshold # threshold = "0.1"
pvalue = opt$pvalue # pvalue = "0.05"
core = opt$core # core = "0.8"
core = as.numeric(core)

output = opt$output

#Create table that will have labels only on interesting points
zipiTableLabels = cbind(zipiTable,zipiTable$Id)
colnames(zipiTableLabels) = c("Id","Pi","Zi","Community","Label")

PiLabels = filter(zipiTableLabels,Pi>=0.62)
ZiLabels = filter(zipiTableLabels,Zi>=2.5)
labeledPoints = rbind(PiLabels,ZiLabels)
unlabeledPoints = subset(zipiTableLabels, !(Id %in% labeledPoints$Id))
unlabeledPoints$Label = ""

zipiTableLabels = rbind(unlabeledPoints,labeledPoints)

#Draw scatterplot
# Write the labels of each quadrant
annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Peripherals","Module hubs"
                   ,"Connectors","Network hubs"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1)) #<- adjust

zipiPlot =  ggplot(zipiTableLabels, aes(x = Pi, y = Zi, label = Label)) +
  geom_point(color = "dark blue", position = "jitter") + geom_hline(yintercept = 2.5, color = "dark green", size = 1) + 
  geom_vline(xintercept = 0.62, color = "orange", size = 1) + 
  labs(x = "Among modules connectivity (Pi)", y = "Within module connectivity (Zi") +
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
  geom_label_repel(max.overlaps = 50) +
  labs(title = paste0("ZiPi plot. # of nodes = ", length(zipiTable$Id),", Core = " ,(core*100),"%, Threshold = ",threshold, ", p-value < ",pvalue)) +
  theme(plot.title = element_text(hjust = 0.5))


gt = ggplot_gtable(ggplot_build(zipiPlot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

ggsave(filename=output,plot=grid::grid.draw(gt))