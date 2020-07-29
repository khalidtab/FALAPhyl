#install.packages("workflow/scripts/vegan_2.5-6.tar", repos = NULL, type="source", INSTALL_opts = '--no-lock')
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="dissimilarity matrix", metavar="Input dissimilarity matrix"),
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

suppressMessages(library("vegan"))
suppressMessages(library("dplyr"))
suppressMessages(library("plotly"))


dissimilarity = opt$input #dissimilarity="/Users/khaled/Desktop/bioinfo_snakemake/data/distance/beta_div/kfp5s2_CPd_GAPd_HNS_LAPd+braycurtis.tsv"
output = opt$output #output="/Users/khaled/Desktop/bioinfo_snakemake/data/plots/svg.svg"  
category = opt$group #category = "ShortCond" 
mapping_file = opt$mapping #mapping_file = "/Users/khaled/Desktop/bioinfo_snakemake/data/map/kfp5s2_CPd_GAPd_HNS_LAPd.txt" 
color = opt$color #color = "Colors" 

dist = read.csv(file=dissimilarity, skip=0, header=T, row.names=1, sep="\t") %>% as.dist(.)
map = read.csv(mapping_file, skip=0, header=T, sep="\t")

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
plotly_coords=merge(coords,map, by.x = "SampleID", by.y = "X.SampleID")
catNum = which(colnames(plotly_coords) == category)
colnames(plotly_coords)[catNum] = "PlotlyCategory"
plotly_coords = plotly_coords %>% dplyr::arrange(PlotlyCategory) # This makes sure that the colors are correctly mapped since the legend is in alphabetical order, and it assumes the entries are in alphabetical order too

myCat = as.vector(t(plotly_coords[catNum]))

colorNum = which(colnames(plotly_coords) == color)
myColor = as.vector(t(plotly_coords[colorNum]))

f1 = list(family = "Arial, sans-serif",size = 15,color = "black")
f2 = list(family = "Arial, sans-serif", size = 15, color = "black")

axis1 = list(title="NMDS 1", titlefont=f1, tickfont=f2, showgrid = T, zeroline = F, gridline = T, linewidth=1, linecolor="#000000")
axis2 = list(title="NMDS 2", titlefont=f1, tickfont=f2, showgrid = T, zeroline = F, gridline = T, linewidth=1, linecolor="#000000")




fig = plot_ly(plotly_coords, x = ~NMDS1, y = ~NMDS2, 
            type = 'scatter', mode = 'markers', 
            name = myCat, size=I(150),
            marker = list(color = myColor, line = list(color = '#000000', width = 1)),            
            alpha=0.75, xaxis=axis1, yaxis=axis2)

for(l in unique(myCat)) {
  myhull = dplyr::filter(plotly_coords,PlotlyCategory == l)
  NMDS1 = as.numeric(as.character(myhull$NMDS1))
  NMDS2 = as.numeric(as.character(myhull$NMDS2))
  hull = chull(x =NMDS1, y = NMDS2)
  hull = c(hull, hull[1])
  myhull = myhull[c(hull),]
  NMDS1 = as.numeric(as.character(myhull$NMDS1))
  NMDS2 = as.numeric(as.character(myhull$NMDS2))
   fig <- fig %>% 
    add_polygons(x=NMDS1,
               y=NMDS2,
               line=list(width=0.5,color="transparent"),
               fillcolor=unique(myhull$Colors), opacity = 0.3, inherit = FALSE, showlegend = FALSE)
}

fig = fig %>% layout(title=paste("stress = ", stress), xaxis=axis1, yaxis = axis2)


myjson = plotly_json(fig,FALSE)
write(myjson,output)
#server$export(fig,file=output)
#server$close()


