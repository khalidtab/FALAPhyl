library("optparse")
suppressMessages(library("dplyr"))
suppressMessages(library("plotly"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Nodes list with Zi-Pi coordinates", metavar="Zi-Pi input"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Name of the output plotly json file", metavar="Plotly json output")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

plotly_coords = read.csv(opt$input, skip=0, header=T, sep="\t")
plotly_coords = plotly_coords %>% mutate(names=case_when(Zi >= 2.5 | Pi >= 0.62 ~ Id))

f1 = list(family = "Arial, sans-serif",size = 15,color = "black")
f2 = list(family = "Arial, sans-serif", size = 15, color = "black")

axis1 = list(title="Among modules connectivity (Pi)", titlefont=f1, tickfont=f2, showgrid = T, zeroline = F, gridline = T, linewidth=1, linecolor="#000000", range = c(-0.1, 0.8))
axis2 = list(title="Within module connectivity (Zi)", titlefont=f1, tickfont=f2, showgrid = T, zeroline = F, gridline = T, linewidth=1, linecolor="#000000", range = c(-2,4))

fig = plot_ly(plotly_coords, x = ~Pi, y = ~Zi, 
            type = 'scatter', mode = 'markers+text', 
             size=I(10),text=~names)#,
        #    marker = list(color = "#000000", line = list(color = '#000000', width = 1)),            
         #   alpha=0.75, xaxis=axis1, yaxis=axis2)


fig = fig %>% layout (xaxis=axis1, yaxis = axis2)

fig = fig %>% add_trace(x = 0.62, showlegend = FALSE, mode = 'lines', size=I(2)) 
fig = fig %>% add_trace(y = 2.5, showlegend = FALSE, mode = 'lines', size=I(2)) 

fig = fig %>% layout (annotations = 
         list(x = 0, y = 0, text = "Peripherals", 
              showarrow = F, xref='paper', yref='paper', 
              xanchor='auto', yanchor='auto', xshift=0, yshift=0,
              font=list(size=12, color="black")))

fig = fig %>% layout (annotations = 
                        list(x = 0, y = 0.85, text = "Module hubs", 
                             showarrow = F, xref='paper', yref='paper', 
                             xanchor='auto', yanchor='auto', xshift=0, yshift=0,
                             font=list(size=12, color="red")))

fig = fig %>% layout (annotations = 
                        list(x = 1.01, y = 0, text = "Connectors", 
                             showarrow = F, xref='paper', yref='paper', 
                             xanchor='auto', yanchor='auto', xshift=0, yshift=0,
                             font=list(size=12, color="blue")))

fig = fig %>% layout (annotations = 
                        list(x = 1.01, y = 0.85, text = "Network hubs", 
                             showarrow = F, xref='paper', yref='paper', 
                             xanchor='auto', yanchor='auto', xshift=0, yshift=0,
                             font=list(size=12, color="green")))

myjson = plotly_json(fig,FALSE)
write(myjson,opt$output)


#plotly::orca(fig,file="~/Desktop/zi-pi/GAP_zi_pi.svg")

