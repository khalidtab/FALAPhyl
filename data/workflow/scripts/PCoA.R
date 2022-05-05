#install.packages("workflow/scripts/vegan_2.5-6.tar", repos = NULL, type="source", INSTALL_opts = '--no-lock')
library("optparse")
suppressMessages(library("vegan"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("cowplot"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Dissimilarity matrix", metavar="Input dissimilarity matrix"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output SVG file", metavar="Output SVG File name"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="The mapping file", metavar="Mapping file"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="The category in the mapping file", metavar="Group name"),
  make_option(c("-x", "--width"), type="character", default=NULL, help="Width of graph", metavar="Width of graph"),
  make_option(c("-y", "--height"), type="character", default=NULL, help="Height of graph", metavar="Height of graph"),
  make_option(c("-c", "--color"), type="character", default=NULL, help="The color column in the mapping file", metavar="Color name")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

dissimilarity = opt$input
output = opt$output
category = opt$group
mapping_file = opt$mapping 
color = opt$color
mywidth = as.numeric(opt$width)
myheight = as.numeric(opt$height)

dissimilarity2 = read_tsv(dissimilarity) %>% suppressMessages(.) %>% as.data.frame(.)
row.names(dissimilarity2) = dissimilarity2[,1]
dissimilarity2[,1] = NULL

PCoA = cmdscale(dissimilarity2,eig=TRUE,k=2) # Only doing 2D graphs 
eigenvaluesSUM = PCoA$eig %>% sum(.)
ProportionsExplainedPCoA1 = (PCoA$eig[1]/eigenvaluesSUM) %>% round(.,4)
ProportionsExplainedPCoA1 = ProportionsExplainedPCoA1*100
ProportionsExplainedPCoA2 = (PCoA$eig[2]/eigenvaluesSUM) %>% round(.,4)
ProportionsExplainedPCoA2 = ProportionsExplainedPCoA2*100


coordinates = PCoA$points %>% as.data.frame(.)
coordinates = cbind(rownames(coordinates),coordinates)
colnames(coordinates) = c("SampleID","PCoA1","PCoA2")

map = suppressMessages(read.csv(mapping_file, skip=0, header=T, sep="\t"))

# Get NMDS coordinates

# Format the NMDS matrix
PCoA_table=merge(coordinates,map, by.x = "SampleID", by.y = "X.SampleID")
catNum = which(colnames(PCoA_table) == category)
colnames(PCoA_table)[catNum] = "Category"
PCoA_table = PCoA_table %>% dplyr::arrange(Category) # This makes sure that the colors are correctly mapped since the legend is in alphabetical order, and it assumes the entries are in alphabetical order too

myCat = as.vector(t(PCoA_table[catNum]))

colorNum = which(colnames(PCoA_table) == color)
myColor = unique(as.vector(t(PCoA_table[colorNum])))

PCoA_table$PCoA1 = as.numeric(as.character(PCoA_table$PCoA1))
PCoA_table$PCoA2 = as.numeric(as.character(PCoA_table$PCoA2))


# Graph with no names
myplot = ggplot(PCoA_table, aes(x = PCoA1, y = PCoA2)) + geom_point(aes(color = Category, size=1)) + 
  labs(title = "Principal coordinates analysis") + theme_bw() + 
  scale_color_manual(values=myColor) + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  # these last 3 lines remove the grid lines
myplot = myplot + labs(x=paste0("PCoA1 - variance explained: ",ProportionsExplainedPCoA1,"%"),y=paste0("PCoA2 - variance explained: ",ProportionsExplainedPCoA2,"%"))

ggsave(filename=output,plot=myplot,width=mywidth,height=myheight)

# Graph with no names  but with probability density plots
myplot = ggplot(PCoA_table, aes(x = PCoA1, y = PCoA2)) + geom_point(aes(color = Category, size=1)) + 
  labs(title = "Principal coordinates analysis") + theme_bw() + 
  scale_color_manual(values=myColor) + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  # these last 3 lines remove the grid lines

xdens = axis_canvas(myplot, axis = "x")+ geom_density(data = PCoA_table, aes(x = PCoA1,fill=Category), alpha = 0.5, size = 0.1)+scale_fill_manual(values=unique(as.vector(t(PCoA_table[colorNum]))))

ydens = axis_canvas(myplot, axis = "y", coord_flip = TRUE)+
  geom_density(data = PCoA_table, aes(x = PCoA2, fill=Category),
               alpha = 0.5, size = 0.1)+coord_flip()+scale_fill_manual(values=unique(as.vector(t(PCoA_table[colorNum]))))
p1 = myplot + xlab(paste0("PCoA1 - variance explained: ",
                          ProportionsExplainedPCoA1,
                      "%"))
p2 = p1 + ylab(paste0("PCoA2 - variance explained: ",
                      ProportionsExplainedPCoA2,
                          "%"))
p3 = insert_xaxis_grob(p2, xdens, grid::unit(.1, "null"), position = "top")
p4 = insert_yaxis_grob(p3, ydens, grid::unit(.1, "null"), position = "right")

ggsave(filename=paste0(output,"noNameswProbDF.svg"),plot=p4,width=mywidth,height=myheight)



# Graph with names
myplot = ggplot(PCoA_table, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(color = Category, size=1)) + 
  labs(title = "Principal coordinates analysis") +
  theme_bw() + 
  scale_color_manual(values=myColor) + 
  geom_text_repel(max.overlaps = 50, label=PCoA_table$SampleID,size=2) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # these last 3 lines remove the grid lines
myplot = myplot + labs(x=paste0("PCoA1 - variance explained: ",ProportionsExplainedPCoA1,"%"),y=paste0("PCoA2 - variance explained: ",ProportionsExplainedPCoA2,"%"))

ggsave(filename=paste0(output,"withnames.svg"),plot=myplot,width=mywidth,height=myheight)

# Graph with names and probability density plots
xdens = axis_canvas(myplot, axis = "x")+ geom_density(data = PCoA_table, aes(x = PCoA1,fill=Category), alpha = 0.5, size = 0.1)+scale_fill_manual(values=unique(as.vector(t(PCoA_table[colorNum]))))

ydens = axis_canvas(myplot, axis = "y", coord_flip = TRUE)+
  geom_density(data = PCoA_table, aes(x = PCoA2, fill=Category),
               alpha = 0.5, size = 0.1)+coord_flip()+scale_fill_manual(values=unique(as.vector(t(PCoA_table[colorNum]))))
p1 = myplot + xlab(paste0("PCoA1 - variance explained: ",
                          ProportionsExplainedPCoA1,
                          "%"))
p2 = p1 + ylab(paste0("PCoA2 - variance explained: ",
                      ProportionsExplainedPCoA1,
                      "%"))
p3 = insert_xaxis_grob(p2, xdens, grid::unit(.1, "null"), position = "top")
p4 = insert_yaxis_grob(p3, ydens, grid::unit(.1, "null"), position = "right")

ggsave(filename=paste0(output,"withnamesprobDF.svg"),plot=p4,width=mywidth,height=myheight)
