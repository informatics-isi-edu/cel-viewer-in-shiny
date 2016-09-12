
library(ggplot2)
library(ggdendro)
library(plotly)
library(RJSONIO)


## make the heatmap data from the json blob

makeHeatmapData_f <- function(jlist) {

#jlist<-fromJSON("../data/HeatmapData.json")

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

  print(vals)
  dimnames(vals) <- list(nsymbols,names)
  hval <- t(vals) # now row=samples, col=symbols
  return(hval)
}

generatePlotlyHeatmap_f <- function(jlist) {

pdf(NULL)
hval <- makeHeatmapData_f(jlist)

#https://plot.ly/r/setting-graph-size/
viewer_m <- list(
  l = 100,
  r = 50,
  b = 100,
  t = 50,
  pad = 4
)

#dendogram data
x <- as.matrix(scale(hval))
dd.col <- as.dendrogram(hclust(dist(x)))
dd.row <- as.dendrogram(hclust(dist(t(x))))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)

# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

# x/y dendograms
px <- ggdend(dx$segments)
py <- ggdend(dy$segments) + coord_flip()

# heatmap
col.ord <- order.dendrogram(dd.col)
row.ord <- order.dendrogram(dd.row)
xx <- scale(hval)[col.ord, row.ord]
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xx)
colnames(df) <- xx_names[[2]]
df$sample <- xx_names[[1]]
mdf <- reshape2::melt(df, id.vars="sample")
p <- ggplot(mdf, aes(x = variable, y = sample)) + 
theme(axis.text.x=element_text(angle=90, hjust=1)) +
     geom_tile(aes(fill = value))


# hide axis ticks and grid lines
eaxis <- list(
  showticklabels = FALSE,
  showgrid = FALSE,
  zeroline = FALSE
)

p_empty <- plotly_empty()
# plot_ly(filename="r-docs/dendrogram") %>%
# note that margin applies to entire plot, so we can
# add it here to make tick labels more readable
#  layout( xaxis = eaxis, yaxis = eaxis)
#pp <- ggplotly(p)


ss <- subplot(px, p_empty, p, py, nrows = 2, margin = 0.02)  %>%
          layout(autosize = F, width = "100%", height = "90%",
                  margin = viewer_m,
  xaxis= list(domain= c(0, 0.80)),
  yaxis= list(domain= c(0.50,1)),
  xaxis2= list(domain= c(0.80, 1)),
  yaxis2= list(domain= c(0.50, 1)),
  xaxis3= list(domain= c(0, 0.80)),
  yaxis3= list(domain= c(0, 0.50)),
  xaxis4= list(domain= c(0.80, 1)),
  yaxis4= list(domain= c(0, 0.50))
                )

  return(ss)
}


## run standalone
##generatePlotlyHeatmap_f(NULL)
