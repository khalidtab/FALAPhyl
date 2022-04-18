#install.packages("workflow/scripts/vegan_2.5-6.tar", repos = NULL, type="source", INSTALL_opts = '--no-lock')
library("optparse")
suppressMessages(library("vegan"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("cowplot"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="dissimilarity matrix", metavar="Input dissimilarity matrix"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output SVG file", metavar="Output SVG File name"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, help="The mapping file", metavar="Mapping file"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="The category in the mapping file", metavar="Group name"),
  make_option(c("-x", "--width"), type="character", default=NULL, help="Width of the SVG", metavar="Width of SVG"),
  make_option(c("-y", "--height"), type="character", default=NULL, help="Height of the SVG", metavar="Height of SVG"),
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


dist = suppressMessages(read.csv(file=dissimilarity, skip=0, header=T, row.names=1, sep="\t")) %>% as.dist(.)
map = suppressMessages(read.csv(mapping_file, skip=0, header=T, sep="\t"))

# Get NMDS coordinates
NMDS = metaMDS(dist)
coords = NMDS$points
stress = NMDS$stress 
stress = round(stress,digits = 3)

# Format the NMDS matrix
firstcol=as.matrix(rownames(coords))
coords=cbind(firstcol,coords)
rownames(coords) = NULL
colnames(coords) = list("SampleID","NMDS1","NMDS2")
NMDS_table=merge(coords,map, by.x = "SampleID", by.y = "X.SampleID")
catNum = which(colnames(NMDS_table) == category)
colnames(NMDS_table)[catNum] = "Category"
NMDS_table = NMDS_table %>% dplyr::arrange(Category) # This makes sure that the colors are correctly mapped since the legend is in alphabetical order, and it assumes the entries are in alphabetical order too

myCat = as.vector(t(NMDS_table[catNum]))

colorNum = which(colnames(NMDS_table) == color)
myColor = unique(as.vector(t(NMDS_table[colorNum])))


NMDS_table$NMDS1 = as.numeric(as.character(NMDS_table$NMDS1))
NMDS_table$NMDS2 = as.numeric(as.character(NMDS_table$NMDS2))

# Graph with no names
myplot = ggplot(NMDS_table, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(color = Category, size=1)) + 
  labs(title = paste("Non-metric Multidimensional Scaling (NMDS), stress = ", stress)) + theme_bw() + 
  scale_color_manual(values=myColor) + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # these last 3 lines remove the grid lines
ggsave(filename=output,plot=myplot,width=mywidth,height=myheight)

# Graph with no names  but with density plots
xdens = axis_canvas(myplot, axis = "x")+ geom_density(data = NMDS_table, aes(x = NMDS1,fill=Category), alpha = 0.5, size = 0.1)+scale_fill_manual(values=unique(as.vector(t(NMDS_table[colorNum]))))

ydens = axis_canvas(myplot, axis = "y", coord_flip = TRUE)+
  geom_density(data = NMDS_table, aes(x = NMDS2, fill=Category),
               alpha = 0.5, size = 0.1)+coord_flip()+scale_fill_manual(values=unique(as.vector(t(NMDS_table[colorNum]))))

p1 = insert_xaxis_grob(myplot, xdens, grid::unit(.1, "null"), position = "top")
p2 = insert_yaxis_grob(p1, ydens, grid::unit(.1, "null"), position = "right")

ggsave(filename=paste0(output,"noNameswithProbDF.svg"),plot=p2,width=mywidth,height=myheight)



# Graph with names
myplot = ggplot(NMDS_table, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(color = Category, size=1)) + 
  labs(title = paste("Non-metric Multidimensional Scaling (NMDS), stress = ", stress)) + theme_bw() + 
  scale_color_manual(values=myColor) + 
  geom_text_repel(max.overlaps = 50, label=NMDS_table$SampleID,size=2) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # these last 3 lines remove the grid lines
ggsave(filename=paste0(output,"withnames.svg"),plot=myplot,width=mywidth,height=myheight)



xdens = axis_canvas(myplot, axis = "x")+ geom_density(data = NMDS_table, aes(x = NMDS1,fill=Category), alpha = 0.5, size = 0.1)+scale_fill_manual(values=unique(as.vector(t(NMDS_table[colorNum]))))

ydens = axis_canvas(myplot, axis = "y", coord_flip = TRUE)+
  geom_density(data = NMDS_table, aes(x = NMDS2, fill=Category),
               alpha = 0.5, size = 0.1)+coord_flip()+scale_fill_manual(values=unique(as.vector(t(NMDS_table[colorNum]))))

p1 = insert_xaxis_grob(myplot, xdens, grid::unit(.1, "null"), position = "top")
p2 = insert_yaxis_grob(p1, ydens, grid::unit(.1, "null"), position = "right")

ggsave(filename=paste0(output,"withnamesProbDF.svg"),plot=p2,width=mywidth,height=myheight)



