
#library(RJSONIO)

leftstr <- function(x, n = 1) substr(x, 1, min(n, nchar(x)))
rightstr <- function(x, n = 1) substr(x, nchar(x) - min(n - 1, nchar(x) - 1), nchar(x))
plot.null <- function()
plot(0, 0, col = "transparent", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
mousecase <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))

#generateMAplotjson_f <- function(ones,twos,serverCFG,dat.sel,dat.top) {
#  MAINLABEL <- paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " ")
#  dir <- gsub(" ", "_", MAINLABEL)
#  dir.create(dir)
generateMAplotjson_f <- function(dat.sel,dat.top) {

#  metaList <- list(type='maplot', config=serverCFG)
  blackX <- dat.sel$A[dat.sel$color == "black"]
  blackY <- dat.sel$M[dat.sel$color == "black"]
  blackSymbol <- dat.sel$symbol[dat.sel$color == "black"]
  blackColor <- dat.sel$color[dat.sel$color == "black"]
  blackPtsList <- list(x=blackX, y=blackY, color=blackColor, symbol=blackSymbol)

  otherX <- dat.sel$A[dat.sel$color != "black"]
  otherY <- dat.sel$M[dat.sel$color != "black"]
  otherSymbol <- dat.sel$symbol[dat.sel$color != "black"]
  otherColor <- dat.sel$color[dat.sel$color != "black"]
  otherPtsList <- list(x=otherX, y=otherY, color=otherColor, symbol=otherSymbol)

  topX <- dat.top$A
  topY <- dat.top$M
  topSymbol <-  dat.top$symbol
  topColor <- dat.top$color
  topPtsList <-list(x=topX, y=topY, color=topColor, symbol=topSymbol)

  dataList <- list(blackPts=blackPtsList)
  if( !is.null(otherX) && length(otherX) != 0) { 
    dataList$otherPts <- otherPtsList
  }
  if(!is.null(topX) && length(topX)!=0) {
    dataList$topPts <- topPtsList
  }
  return( list(data=dataList) )
##  jsonList <- list(meta=metaList, data=dataList)
##  maplotfn <- paste0(dir, "/", "MAplotData.json")
##  write(toJSON(jsonList), maplotfn, append=FALSE)
##  return(jsonList)
}

#generateHeatmapjson_f <- function(ones,twos,serverCFG,dat.top,dat.heat) {
#  MAINLABEL <- paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " ")
#  dir <- gsub(" ", "_", MAINLABEL)
#  dir.create(dir)
generateHeatmapjson_f <- function(topsymbol,dat.heat) {

  SYMBOL <- topsymbol
  PROBESET <- rownames(dat.heat) 

#  metaList <- list(type='heatmap', config=serverCFG)
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
  return (list( data=heatmapList) )
##  jsonList <- list(meta=metaList, data=heatmapList)
##  heatmapfn <- paste0(dir, "/", "HeatmapData.json")
##  write(toJSON(jsonList),heatmapfn, append=FALSE) 
##  return(jsonList)
}


processConfig_f <- function(query) {

## set default, inverted is initalized to false
  serverCFG <- list( sel=c("E10.5_Mnd_D","E10.5_Mnd_P"),
                     comp= c("place"),
                     url= NULL,
                     summary =c("Z"),
                     gene1 = "dumdum",
                     log = T,
                     fc=2,
                     fdr=0.01,
                     inverted = "normal",
                     max=20,
                     ages=10.5
                   )

  if (!is.null(query[['url']])) {
    serverCFG$url <- query[['url']]

  } 
  if (!is.null(query[['sel']]))  {
# stert with "E10.5_Mnd_D E10.5_Mnd_P"
# into  "E10.5_Mnd_D" "E10.5_Mnd_P"
      s <- query[['sel']]
      s <- gsub("\"", "", s)
      ss <- strsplit(s," ")
      ss <- ss[[1]]
      t <- character()
      for(s in ss) t <- c(t, s)
      serverCFG$sel <- t
      tsel <- substr(t, 2, nchar(t))
      tsel <- gsub("_", "", tsel)
      tsels <- character()
      for (t in tsel) tsels <- c(tsels, paste0(t, 1:3))
      serverCFG$sels <- tsels
      ages <- unique(substr(tsels, 1, 4))
      serverCFG$ages <-ages 
    } 

    if (!is.null(query[['comp']])) { 
      s <- query[['comp']]
      s <- gsub("\"", "", s)
      serverCFG$comp <- c(s)
    }

    ## do the selection processing

    sels <- serverCFG$sels

    if (serverCFG$comp == "place"){ 
      serverCFG$target.col <- "placeD"
      ones <- grep("P", sels, value = T)
      twos <- grep("D", sels, value = T)
      ones.reduced <- gsub("P", "", ones)
      twos.reduced <- gsub("D", "", twos)
      serverCFG$ones <- ones[ones.reduced %in% twos.reduced]
      serverCFG$twos <- twos[twos.reduced %in% ones.reduced]
      invert <- ifelse(serverCFG$inverted == "inverted", T, F)
      serverCFG$one <- ifelse(invert, "distal", "proximal")
      serverCFG$two <- ifelse(invert, "proximal", "distal")
    }
    if (serverCFG$comp == "bone"){
      serverCFG$target.col <- "boneMnd"
      ones <- grep("Max", sels, value = T)
      twos <- grep("Mnd", sels, value = T)
      ones.reduced <- gsub("Max", "", ones)
      twos.reduced <- gsub("Mnd", "", twos)
      serverCFG$ones <- ones[ones.reduced %in% twos.reduced]
      serverCFG$twos <- twos[twos.reduced %in% ones.reduced]
      invert <- ifelse(serverCFG$inverted == "inverted", T, F)
      serverCFG$one <- ifelse(invert, "maxilla", "mandible")
      serverCFG$two <- ifelse(invert, "mandible", "maxilla")
    }

    if (serverCFG$comp == "age"){
      ages <- serverCFG$ages
      if (length(ages)) age.1 <- min(ages) else age.1 <- 10.5
      age.2 <- max(ages)
      serverCFG$target.col <- paste0("age", age.2)
      ones <- grep(age.1, sels, value = T)
      twos <- grep(age.2, sels, value = T)
      ones.reduced <- gsub(age.1, "", ones)
      twos.reduced <- gsub(age.2, "", twos)
      serverCFG$ones <- ones[ones.reduced %in% twos.reduced]
      serverCFG$twos <- twos[twos.reduced %in% ones.reduced]
      invert <- ifelse(serverCFG$inverted == "inverted", T, F)
      serverCFG$one <- ifelse(invert, paste0("E", age.2), paste0("E", age.1))
      serverCFG$two <- ifelse(invert, paste0("E", age.1), paste0("E", age.2))
    }
    return(serverCFG)
}

processDatSel_f <- function(serverCFG, dat) {

## get preset 
  invert <- ifelse(serverCFG$inverted == "inverted", T, F)
  ones <- serverCFG$ones
  twos <- serverCFG$twos

  genes.tab <- xtabs(~unlist(dat$symbol))
  single.genes <- names(genes.tab)[genes.tab == 1]

#### 
#### HERE is when the data is processed from full raw set to raw selected set
####
  dat.sel <- dat[, c(ones, twos, "symbol")]

  if (serverCFG$summary == "AVE"){
    dat.sel <- ddply(dat.sel, .(symbol), numcolwise(mean))
    dat.sel <- dat.sel[dat.sel$symbol != "NA", ]
    rownames(dat.sel) <- dat$probeset[match(dat.sel$symbol, dat$symbol)]}
  if (serverCFG$summary == "A"){
    dat.sel$A <- rowMeans(dat.sel[, c(ones,twos)])
    dat.sel$P <- rownames(dat.sel)
    dat.single <- dat.sel[dat.sel$symbol %in% single.genes, ]
    dat.multiple <- dat.sel[!dat.sel$symbol %in% single.genes, ]
    dat.agg <- aggregate(A ~ symbol, dat.multiple, max)
    dat.multiple <- merge(dat.agg, dat.multiple, sort = F)
    rownames(dat.multiple) <- dat.multiple$P
    dat.sel <- rbind(dat.single, dat.multiple)
    dat.sel <- dat.sel[, colnames(dat.sel)[!colnames(dat.sel) %in% c("A", "P")]]}
  if (serverCFG$summary == "M"){
    dat.sel$MM <- abs(rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones]))
    dat.sel$P <- rownames(dat.sel)
    dat.single <- dat.sel[dat.sel$symbol %in% single.genes, ]
    dat.multiple <- dat.sel[!dat.sel$symbol %in% single.genes, ]
    dat.agg <- aggregate(MM ~ symbol, dat.multiple, max)
    dat.multiple <- merge(dat.agg, dat.multiple, sort = F)
    rownames(dat.multiple) <- dat.multiple$P
    dat.sel <- rbind(dat.single, dat.multiple)
    dat.sel <- dat.sel[, colnames(dat.sel)[!colnames(dat.sel) %in% c("MM", "P")]]}
    return(dat.sel)
}

makeTop_f <- function(dat,serverCFG,design) {

## get preset 
    invert <- ifelse(serverCFG$inverted == "inverted", T, F)
    ones <- serverCFG$ones
    twos <- serverCFG$twos
    one <- serverCFG$one
    two <- serverCFG$two
    target.col <- serverCFG$target.col
    sels <- serverCFG$sels
    gene1 <- serverCFG$gene1
    log <- serverCFG$log
    fc <- serverCFG$fc
    fdr <- serverCFG$fdr
    max <- serverCFG$max

    if (!nrow(design) | !target.col %in% colnames(design) | !all(sels %in% c(ones, twos))) {
      return( list() )
    }

    dat.sel <- processDatSel_f(serverCFG, dat)

    fit <- lmFit(dat.sel[, colnames(dat.sel)[colnames(dat.sel) != "symbol"]], design)
    efit <- eBayes(fit)

    dat.sel$A <- rowMeans(dat.sel[, colnames(dat.sel)[colnames(dat.sel) != "symbol"]])
    dat.sel$M <- rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones])
    if (invert) dat.sel$M <- -dat.sel$M
    if (!is.null(dat.sel$symbol))
      dat.sel$symbol[is.na(dat.sel$symbol)] <-  rownames(dat.sel)[is.na(dat.sel$symbol)]
    dat.sel$color <- "black"

## initialize that the gene1 to be empty
    if(gene1 != "") {
      gene1 <- c(strsplit(gene1," "))
      gene1 <- gene1[[1]]
      if( length(gene1) > 0 ) {
        for (g in gene1) {
          dat.sel$color[dat.sel$symbol == mousecase(g) | rownames(dat.sel) == mousecase(g)] <- "blue"
        }
      }
    }

    lfc <- ifelse(log, abs(as.numeric(fc)), log2(abs(as.numeric(fc))))
    top <- topTable(efit, coef = which(colnames(design) == target.col), 
      lfc = lfc, p.value = as.numeric(fdr), num = Inf)
    top <- merge(top, dat[, c("probeset", "symbol", "desc")], by.x = 0, by.y = "probeset")
    top.colored <- top[rownames(top) %in% dat.sel$probeset[dat.sel$color != "black"], ]

# top.colored seems to have nrow of 0 even though ncol has something
# ???
    if (length(top.colored)){
      top.black <- top[1:min(nrow(top), as.numeric(max)), ]
      top <- merge(top.colored, top.black, all = T)
    }
    if (nrow(top)){
      if (invert) top$logFC <- -top$logFC      
      top <- top[order(top$logFC, decreasing = T), ]
      top$FC <- 2 ^ top$logFC
      top$FC[top$logFC < 0] <- -1 / top$FC[top$logFC < 0]
      top$AveExpr <- round(top$AveExpr, 3)
      top$logFC <- round(top$logFC, 2)
      top$FC <- round(top$FC, 2)
      top$adj.P.Val <- format(top$adj.P.Val, scientific = T, digits = 2)
      top$desc <- format(top$desc, justify = "left")
      top$symbol <- format(top$symbol, justify = "left")
      top <- top[, c("Row.names", "symbol", "AveExpr", "logFC", "FC", "adj.P.Val", "desc")]
      colnames(top) <- c("Probeset", "Gene", "Ave.Expr", "Log2FC", "Fold.Change", "Adj.P.Val", "Description")
    } else top <- data.frame(NULL)
    dat.top <- dat.sel[rownames(dat.sel) %in% top$Probeset, ]

    tlist <- list(tt=top, dd=dat.top, ss=dat.sel)
    return(tlist)
}

makeDesignControl_f <- function(dat, serverCFG) {

    ones <- serverCFG$ones
    twos <- serverCFG$twos
    age.1 <- min(serverCFG$ages)
    target.col <- serverCFG$target.col

    age <- factor(rep(c(10.5, 11.5, 12.5, 13.5, 14.5), each = 12))
    bone <- factor(rep(rep(c("Max", "Mnd"), each = 6), 5))
    place <- relevel(factor(rep(rep(c("D", "P"), each = 3), 10)), "P")
    full.design <- model.matrix(~age + bone + place)
    rownames(full.design) <- paste0(age, bone, place, 1:3)

    young.col <- paste0("age", age.1)

    full.des <- full.design[c(ones, twos), ]
    constant <- apply(full.des, 2, function(x) length(unique(x)) == 1)
    constant[1] <- F
    constant[names(constant) == young.col] <- T
    design.cols <- names(constant[constant == F])
    design <- data.frame(lapply(design.cols, function(x){col <- list(full.des[, x])}))
    names(design) <- design.cols
    control.cols <- design.cols[design.cols != target.col]
    control <- data.frame(lapply(control.cols, function(x){col <- list(full.des[, x])}))
    names(control) <- control.cols
    plist <- list(dd=design, cc=control)

    return(plist)
}

makeHeat_f <- function(dat.top, serverCFG, control ) {

  sels <- serverCFG$sels
  heatadjust <- serverCFG$heatadjust
  heatscale <- serverCFG$heatscale

  dat.heat <- as.matrix(dat.top[, sels])
  dat.heat <- dat.heat - mean(dat.heat)
  if (heatadjust){
    for (col in colnames(control)[colnames(control) != "(Intercept)"]){
      for (v in unique(control[, col])){
        mm <- mean(dat.heat[, rownames(control)[control[, col] == v]])
              dat.heat[, rownames(control)[control[, col] == v]] <-
                dat.heat[, rownames(control)[control[, col] == v]] - mm}}}
  if (heatscale %in% c("MC", "Z")) dat.heat <- sweep(dat.heat, 1, rowMeans(dat.heat))
  if (heatscale == "Z") dat.heat <- sweep(dat.heat, 1, apply(dat.heat, 1, sd), "/")

  return(list(hh=dat.heat))
}


