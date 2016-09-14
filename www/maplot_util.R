library(RJSONIO)
library(plotly)


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


makeMAplotData_f <- function(jlist) {

## DEBUG-IT
jlist<-fromJSON("../data/MAplotData.json")
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
  i<-0
  sz<-length(topY)
  while(i<=sz) {
      if(topY[[i]] >= 0) {
          posX <- c(posX, topX[[i]])
          posY <- c(posY, topY[[i]])
          posColor <- c(posColor, topColor[[i]])
          posSymbol <- c(posSymbol, topSymbol[[i]])
      } else {
          negX <- c(negX, negX[[i]])
          negY <- c(negY, negY[[i]])
          negColor <- c(negColor, negColor[[i]])
          negSymbol <- c(negSymbol, negSymbol[[i]])
      }
      i=i+1
  }
  posPts <- list( x=posX, y=posY, symbol=posSymbol, color=posColor)
  negPts <- list( x=negX, y=negY, symbol=negSymbol, color=negColor)
  newdata <- list( blackPts = blackPts, posPts=posPts, negPts=negPts)
  return(newdata)
}


## format is a list
#  {
#  "blackPts": { "x": [ 10.317,..], "y": [ ...], "symbol":[..], "color":[...] },
#  "posPts": { "x": [ 10.317,..], "y": [ ...], "symbol":[..], "color":[...] },
#  "negPts": { "x": [ 10.317,..], "y": [ ...], "symbol":[..], "color":[...]}
#  }
generatePlotlyMAplot_f <- function(jlist) {

pdf(NULL)
hval <- makeMAplotData_f(jlist)

## figure out the range of y from blackPts+posPts+negPts' 
blackPts <- hval$blackPts
posPts <- hval$posPts
negPts <- hval$negPts
 
allY <- c( blackPts$y, posPts$y, negPts$y)
lim <- max(abs(range(allY)))
print(lim)
}

## DEBUG-IT, run standalone
p <- generatePlotlyMAplot_f(NULL)
print(p)
                                 
