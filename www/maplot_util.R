library(plotly)

##DEBUG-IT
##library(RJSONIO)


## make the maplot data from the json blob
## input json format
## {
## "meta": {
##      "title": "10.5MndP 10.5MndD",
##      "config": { // from configIt.R }
##         },
##  "data" : {
##      "blackPts": { "x":[ 10.317,..],"y": [ ...],"symbol":[..],"color":[...] },
##      "otherPts": { "x":[ 10.317,..],"y": [ ...],"symbol":[..],"color":[...] },
##      "topPts": { "x":[ 10.317,..],"y": [ ...],"symbol":[..],"color":[...] }
##           }
## }

grColor <- c('rgb(0,139,0)', 'rgb(139,0,0)')
rgColor <- c('rgb(139,0,0)', 'rgb(0,139,0)')

makeMAplotData_f <- function(jlist) {

##DEBUG-IT
##jlist<-fromJSON("MAplotData.json")
  data <- jlist$data
  blackPts <- data$blackPts
  otherPts <- data$otherPts
  topPts <- data$topPts
## split topPts into positive and negative
  posX <- c()
  posY <- c()
  posColor <- c()
  posSymbol <- c()
  negX <- c()
  negY <- c()
  negColor <- c()
  negSymbol <- c()
  topX <- topPts$x
  topY <- topPts$y
  topColor <- topPts$color
  topSymbol <- topPts$symbol
  i<-1
  sz<-length(topY)
  while(i<=sz) {
      if(topY[[i]] >= 0) {
          posX <- c(posX, topX[[i]])
          posY <- c(posY, topY[[i]])
          posColor <- c(posColor, topColor[[i]])
          posSymbol <- c(posSymbol, topSymbol[[i]])
      } else {
          negX <- c(negX, topX[[i]])
          negY <- c(negY, topY[[i]])
          negColor <- c(negColor, topColor[[i]])
          negSymbol <- c(negSymbol, topSymbol[[i]])
      }
      i=i+1
  }
  posPts <- list( x=posX, y=posY, symbol=posSymbol, color=posColor)
  negPts <- list( x=negX, y=negY, symbol=negSymbol, color=negColor)
  newdata <- list( blackPts = blackPts, posPts=posPts, negPts=negPts, otherPts=otherPts)
  return(newdata)
}

## format is a list
#  {
#  "blackPts": { "x": [ 10.317,..], "y": [ ...], "symbol":[..], "color":[...] },
#  "posPts": { "x": [ 10.317,..], "y": [ ...], "symbol":[..], "color":[...] },
#  "negPts": { "x": [ 10.317,..], "y": [ ...], "symbol":[..], "color":[...]}
#  "otherPts": { "x": [ 10.317,..], "y": [ ...], "symbol":[..], "color":[...]}
#  }
#
#  GenesTrace=3; // 0=blackPts, 1=line(2), 2=line(-2) 3=posPts, 4=negPts, 5=special
#
generatePlotlyMAplot_f <- function(jlist,my.xlab,my.ylab,my.title,my.lfc) {

viewer_m <- list(
  l = 100,
  r = 50,
  b = 100,
  t = 50,
  pad = 4
)

pdf(NULL)
hval <- makeMAplotData_f(jlist)

## figure out the range of y from just blackPts
blackPts <- hval$blackPts
posPts <- hval$posPts
negPts <- hval$negPts
otherPts <- hval$otherPts
 
## with alittle cushion
allY <- c( blackPts$y, posPts$y, negPts$y, otherPts$y)
ylim <- ceiling(max(abs(range(allY))) * 1.2)
yrange_max <- ylim
yrange_min <- ylim * (-1)

allX <- c( blackPts$x)
cntBlacks <- ceiling(length(allX) * 0.1)
xrange_max <- ceiling(max(range(allX)))+1
xrange_min <- floor(min(range(allX)))-1
##

markerlist <- list(size=6, color='#2e2e2e', symbol=100)

blacks <- c()
blacks$x <- head(blackPts$x,cntBlacks)
blacks$y <- head(blackPts$y,cntBlacks)
blacks$symbol <- head(blackPts$symbol,cntBlacks)

p <- plot_ly(data=blacks, type='scatter', x=x, y=y, 
                  textposition="top right",
                  text=symbol,
                  hoverinfo="x+y+text",
                  showlegend = FALSE,
                  marker=markerlist,  mode="markers") %>%
  layout( hovermode = 'closest',
          margin = viewer_m, 
          title = my.title,
          xaxis = list( title=my.xlab, range= list(xrange_min, xrange_max)),
          yaxis = list( title=my.ylab, range= list(yrange_min, yrange_max)))

xlist <- xrange_min:xrange_max
p <- addMAplotLineTrace_f(p, xlist=xlist, yval=my.lfc, 'grey')
p <- addMAplotLineTrace_f(p, xlist=xlist, yval=(-1 *my.lfc), 'grey')

if( ! is.null(otherPts) ) {
  p <- addMAplotDataTrace_f(p, posPts, "top center", '#008b00')
  p <- addMAplotDataTrace_f(p, negPts, "top center", '#8b0000')
  p <- addMAplotDataTrace_f(p, otherPts, "top center", '#00FFFF')
  } else {
    p <- addMAplotDataTrace_f(p, posPts, "top center", '#008b00')
    p <- addMAplotDataTrace_f(p, negPts, "bottom center", '#8b0000')
}

return (p)
}

# "Pts": { "x": [ 10.317,..], "y": [ ...], "symbol":[..], "color":[...]}
# add a trace to a scatter plot
addMAplotDataTrace_f <- function(p, nPts, nPos, nColor) {
  markerlist <- list(size=8,color=nColor)
  q <- p %>% add_trace(data=nPts, type='scatter', x=x, y=y, 
                  textposition=nPos,
                  showlegend = FALSE,
                  text=symbol,
                  hoverinfo="x+y+text",
                  marker=markerlist,  mode="text+markers")
  return(q)
}

addMAplotLineTrace_f <- function(p, xlist, yval, nColor) {
  n=length(xlist)
  ylist <- rep(yval:yval,each=n)
  q <- p %>% add_trace(x = xlist, y = ylist,
        opacity = 0.4,
        line = list(dash="dashed", color=nColor),
        showlegend = FALSE,
        mode="lines")
  return(q)
}

##DEBUG-IT, run standalone
##p <- generatePlotlyMAplot_f(jlist=NULL,my.xlab="X",my.ylab="Y", my.title="ABC")
##print(p)
                                 
