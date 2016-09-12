
library(RJSONIO)

## load in the raw gene expression data
## 'dat' is the dataframe's name
loadGeneDataViaURL_f <- function(inputURL) {
  repeat{
      con <- url(inputURL)
      t <- try(load(con))
      close(con)
      if (!inherits(t, "try-error")) break}
}

generateMAPLOTjson_f <- function(ones,twos,inputCONFIG,dat.sel,dat.top) {
  MAINLABEL <- paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " ")
  dir <- gsub(" ", "_", MAINLABEL)
  dir.create(dir)

  metaList <- list(type='maplot', config=inputCONFIG)
  blackX <- dat.sel$A[dat.sel$color == "black"]
  blackY <- dat.sel$M[dat.sel$color == "black"]
  blackSymbol <- dat.sel$symbol[dat.sel$color == "black"]
  blackColor <- dat.sel$color[dat.sel$color == "black"]
  blackPtsList <- list(x=blackX, y=blackY, symbol=blackSymbol, color=blackColor)

  otherX <- dat.sel$A[dat.sel$color != "black"]
  otherY <- dat.sel$M[dat.sel$color != "black"]
  otherSymbol <- dat.sel$symbol[dat.sel$color != "black"]
  otherColor <- dat.sel$color[dat.sel$color != "black"]
  otherPtsList <- list(x=otherX, y=otherY, color=otherColor, symbol=otherSymbol)

  topX <- dat.top$A
  topY <- dat.top$M
  topSymbol <-  dat.top$symbol
  otherColor <- dat.top$color
  topPtsList <-list(x=topX, y=topY, symbol=topSymbol)

  dataList <- list(blackPts=blackPtsList)
  if( !is.null(otherX) && length(otherX) != 0) { 
    dataList$otherPts <- otherPtsList
  }
  if(!is.null(topX) && length(topX)!=0) {
    dataList$topPts <- topPtsList
  }
  jsonList <- list(meta=metaList, data=dataList)
  maplotfn <- paste0(dir, "/", "MAplotData.json")
##  write(toJSON(jsonList), maplotfn, append=FALSE)
  return(jsonList)
}

generateHEATMAPjson_f <- function(ones,twos,inputCONFIG,dat.top,dat.heat) {
  MAINLABEL <- paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " ")
  dir <- gsub(" ", "_", MAINLABEL)
  dir.create(dir)

  SYMBOL <- dat.top$symbol
  PROBESET <- rownames(dat.heat) 

  metaList <- list(type='heatmap', config=inputCONFIG)
  heatmapList <- list(symbol=SYMBOL, probeset=PROBESET) 

  sampleList.names <- colnames(dat.heat)
  sampleList <- vector("list", length(sampleList.names))
  N <- ncol(dat.heat)
  for(i in as.numeric(1:N)) {
     c <- sampleList.names[[i]]
     cc<-dat.heat[,c]
     ccc <- as.vector(cc)
     l <- list( name=c, data=ccc)
     sampleList[[i]] <- l
  }

  heatmapList$samples <- sampleList
  jsonList <- list(meta=metaList, data=heatmapList)
  heatmapfn <- paste0(dir, "/", "HeatmapData.json")
#  write(toJSON(jsonList),heatmapfn, append=FALSE) 
  return(jsonList)
}
