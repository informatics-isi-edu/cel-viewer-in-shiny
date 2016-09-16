library(ggdendro)
library(reshape2)
#library(RJSONIO)


## make the heatmap data from the json blob
##  json blob format
## 
##  {
##  "meta": {
##      "type": "heatmap"
##      "config" : { from configIt.R }
##  }
##  "data": {
##      "symbol": [...],
##      "probeset": [ ...],
##      "samples": [{ "name":"10.5M.", "data":[...]},
##                  { "name":"... } ]
##  }
##  } 
makeHeatmapData_f <- function(jlist) {

## DEBUG-IT
##jlist<-fromJSON("../data/HeatmapData.json")

  data <- jlist$data
  symbol <- data$symbol ## genes
  sz <- length(symbol)
  nsymbols <- c()
  i<-1
  while(i<=sz) {
      n <- paste0(symbol[[i]],'_',i)
      nsymbols<-c(nsymbols,n)
      i=i+1
  }
  samples <- data$samples
  names<-sapply(samples, function(x) x[[1]])
  vals<-sapply(samples, function(x) x[[2]])

  dimnames(vals) <- list(nsymbols,names)
  hval <- t(vals) # now row=samples, col=symbols
  return(hval)
}

generatePlotlyHeatmap_f <- function(jlist, distfun.row, distfun.col) {

pdf(NULL)
hval <- makeHeatmapData_f(jlist)

#https://plot.ly/r/setting-graph-size/
viewer_m <- list(
  l = 100,
  r = 10,
  b = 100,
  t = 20,
  pad = 2
)

#dendogram data
## row=genes, col=samples
#x <- as.matrix(scale(hval))
# x is 6x22
#dd.col <- as.dendrogram(hclust( dist(x) ))
#dd.col, with 2 branches and 6 members, at height 8.781363
#dd.row <- as.dendrogram(hclust( cor.dist( t(x), abs=F )))
#dd.row, with 2 branches and 22 members, at height 0.05241376
x <- as.matrix(hval)
dd.col <- as.dendrogram(hclust( distfun.col(x) ))
dd.row <- as.dendrogram(hclust( distfun.row( t(x) )))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)

# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.border = element_blank(), 
          panel.grid = element_blank())
}

# x/y dendograms
px <- ggdend(dx$segments)
py <- ggdend(dy$segments) + coord_flip()

# heatmap
col.ord <- order.dendrogram(dd.col)
row.ord <- order.dendrogram(dd.row)
## DEBUG-IT
##browser()

xx <- hval[col.ord, row.ord]
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xx)
colnames(df) <- xx_names[[2]]

df$sample <- xx_names[[1]]
mdf <- reshape2::melt(df, id.vars="sample")
## change from, [1] "sample"   "variable" "value"   
## to, [1] "sample" "gene"   "value" 
names(mdf)<-c("sample","gene","value")
p <- ggplot(mdf, aes(x = gene, y = sample, fill= value)) +
     labs(x = "", y = "", fill = "") +
     geom_tile() +
     theme_bw() +
     theme(axis.text.x=element_text(angle=45, hjust=1))+
     theme(panel.border=element_blank()) +
     scale_fill_gradient2(low = '#FF0000', mid='#000000', high = '#00FF00')


# hide axis ticks and grid lines
eaxis <- list( showgrid = FALSE, zeroline = FALSE )

p_empty <- plotly_empty()
# plot_ly(filename="r-docs/dendrogram") %>%
# note that margin applies to entire plot, so we can
# add it here to make tick labels more readable
#  layout( xaxis = eaxis, yaxis = eaxis)
#pp <- ggplotly(p)


ss <- subplot(px, p_empty, p, py, nrows = 2, margin=0, 
      shareX=TRUE, shareY=TRUE, widths=c(0.8,0.2), heights=c(0.3,0.7))  %>%
      layout( margin = viewer_m)
  return(ss)
}

## DEBUG-IT, run standalone
##p <- generatePlotlyHeatmap_f(NULL)
##print(p)
