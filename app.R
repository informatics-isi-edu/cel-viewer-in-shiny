library(bioDist)
library(gplots)
library(gtools)
library(limma)
library(DT)
library(png)
library(plotrix)
library(plyr)

library(plotly)

library(shiny)
options(shiny.trace = TRUE)
options(shiny.fullstacktrace = TRUE)
options(shiny.error=traceback)

library(magrittr)
library(shinyjs)
library(logging)

basicConfig()

options(shiny.error = function() { 
    logging::logerror(sys.calls() %>% as.character %>% paste(collapse = ", ")) })

inputCONFIG <- list( sel=c("E10.5_Mnd_D","E10.5_Mnd_P"), comp= c("place"))

ui <- fluidPage(
useShinyjs(),
titlePanel("Gene Expression Comparison App"),
sidebarLayout(
  div(class = "col-sm-4",
  tags$form(class = "well",
    textInput("url", "URL", "https://dev.facebase.org/~mei/data/urlexprs.R"),
#    textInput("url", "URL", "data/exprs.R"),
    div(style="padding-bottom:3px", "Chai, Yang; Parada, Carolina; Grimaldi, Alexandre; Ho, Thach-Vu; Harunaga, Jill; Samuels, Bridget; Johnson, Daniel."),
    div(style="padding-bottom:3px", "Distal/proximal mandible/maxilla at embryonic day 10.5/11.5/12.5/13.5/14.5, 3 replicates/condition."),
    div(style="padding-bottom:3px", "Affymetrix Mouse Genome 430 2.0 Array."),
    div(style="padding-bottom:7px", "Mouse drawing by Darrin Lunde and Nguyen Truong Son."),
    conditionalPanel("!output.datloaded",
        helpText("Data is loading (may take a minute) ...")),
    conditionalPanel("output.datloaded", 
        helpText("Data is loaded."))),
  tags$form(class = "well",
      div(tags$label("Highlight one or more genes, ")),
      bootstrapPage(
      div(style="display:inline-block; width:100%", textInput("gene1", NULL, "name")),
      div(style="padding-bottom:3px", "(eg. Tnnt2 Myl3 Rgs5)"))
      ),

  tags$form(class = "well",
   if (inputCONFIG$comp == "place"){ 
     div(tags$label("Select direction for comparing place,"))
   } else if (inputCONFIG$comp == "bone") {
     div(tags$label("Selection direction for comparing bone,"))
   } else if (inputCONFIG$comp == "age") {
     div(tags$label("Selection direction for comparing age,"))
   },
   bootstrapPage(
   if (inputCONFIG$comp == "place"){ 
       radioButtons("invert_place", br(), 
           c("proximal down, distal up" = "normal",
             "distal down, proximal up" = "inverted"), "normal", inline = T)
   } else if (inputCONFIG$comp == "bone") {
       radioButtons("invert_bone", br(), 
           c("mandible down, maxilla up" = "normal",
             "maxilla down, mandible up" = "inverted"), "normal", inline = T)
   } else if (inputCONFIG$comp == "age") {
       radioButtons("invert_age", br(),
           c("earliest down, latest up" = "normal",
             "latest down, earliest up" = "inverted"), "normal", inline = T)
   }),
   conditionalPanel("output.numsel > 6",
       checkboxInput("heatadjust", "For heatmap, factor out comparisons other than the one selected", value = T))),
  tags$form(class= "well",
      bootstrapPage(
      div(style="display:inline-block; width:30%", textInput("fc", "Fold change cut-off", 2)),
        div(style="display:inline-block", checkboxInput("log", "log2", T))),    
    textInput("fdr", "False discovery rate", 0.01, width = "30%"),
    textInput("max", "Max. DE genes/probesets", Inf, width = "30%"),
    radioButtons("summary", "For each gene, show",
        choices = c("average probeset" = "AVE",  "most highly expressed probeset" = "A",
          "most differentially expressed probeset" = "M", "all probesets" = "Z"),
          selected = "AVE", inline = T),
    radioButtons("heatcol", "For heatmap, high values are",
        choices = c("green" = "rg", "red" = "gr"), selected = "rg", inline = T),
    radioButtons("distfun", "For heatmap, cluster genes using",
        choices = c("absolute correlation" = "AC", "correlation" = "C", "Euclidean distance" = "E"),
          selected = "AC", inline = T),
    radioButtons("heatscale", "For heatmap, scale genes by",
        choices = c("mean-centering" = "MC", "Z-score" = "Z", "none" = "N"), selected = "MC", inline = T))),

  mainPanel(width = 8,

    div(style="border:2px solid orange", plotOutput("ma.plot", height = "500px")),
    div(style="border:2px solid blue", plotlyOutput("ma.plotly", height = "500px")),

    div(style="border:2px solid green", plotOutput("heatmap", height = "450px")),
    div(style="border:2px solid red", plotlyOutput("heatmapPlotly", height = "450px")),

    downloadButton("download.table", "Download table"),  br(), br(),
    dataTableOutput("table", width = "90%"))
) ## sidebarLayout
)

server <- function(input, output, session){

printLogJs <- function(x, ...) { 
                               logjs(x) 
                               T
                               }
addHandler(printLogJs)

## process args that were passed in via url

##https://dev.facebase.org/shiny/apps/url/?url=hello-world
##
##https://dev.facebase.org/shiny/apps/usc/?url=https://dev.facebase.org/~mei/data/urlexprs.R
  source("www/old_util.R")
  source("www/util.R")
  source("www/heatmap_util.R")
  source("www/maplot_util.R")

  load("data/jaws3.R")

#  load("data/exprs.R")
#  repeat{
#      con <- url(input$url)
#      con <- input$url
#      if( con != "data/exprs.R") {
#        load("data/exprs.R")
#      }
#  if (!inherits(t, "try-error")) break}
#  myurl <- input$url
#  loginfo("here..")
#  cat("here message\n", file=stderr())
#  load("data/exprs.R")
#  query <- parseQueryString(session$clientData$url_search)
#  if (!is.null(query[['text']])) {
#    q <- query[['text']]
#    updateTextInput(session, "url", value = q)
#  }

  con <- url("https://dev.facebase.org/~mei/data/urlexprs.R")
#  con <- url(input$url)
  load(con)
  

  pdf(NULL)

  age <- factor(rep(c(10.5, 11.5, 12.5, 13.5, 14.5), each = 12))
  bone <- factor(rep(rep(c("Max", "Mnd"), each = 6), 5))
  place <- relevel(factor(rep(rep(c("D", "P"), each = 3), 10)), "P")
  full.design <- model.matrix(~age + bone + place)
  rownames(full.design) <- paste0(age, bone, place, 1:3)

  rownames(dat) <- dat[, "probeset"]
  genes.tab <- xtabs(~unlist(dat$symbol))
  single.genes <- names(genes.tab)[genes.tab == 1]
  output$datloaded <- reactive(is.numeric(nrow(dat)))
  outputOptions(output, "datloaded", suspendWhenHidden = FALSE)
  
  output$table <- renderDataTable({
    sel <- inputCONFIG$sel
    sel <- substr(sel, 2, nchar(sel))
    sel <- gsub("_", "", sel)
    sels <- character()
    for (s in sel) sels <- c(sels, paste0(s, 1:3))
    if (inputCONFIG$comp == "place"){ 
      target.col <- "placeD"
      ones <- grep("P", sels, value = T)
      twos <- grep("D", sels, value = T)
      ones.reduced <- gsub("P", "", ones)
      twos.reduced <- gsub("D", "", twos)
      ones <- ones[ones.reduced %in% twos.reduced]
      twos <- twos[twos.reduced %in% ones.reduced]
      invert <- ifelse(input$invert_place == "inverted", T, F)
      one <- ifelse(invert, "distal", "proximal")
      two <- ifelse(invert, "proximal", "distal")}
    if (inputCONFIG$comp == "bone"){
      target.col <- "boneMnd"
      ones <- grep("Max", sels, value = T)
      twos <- grep("Mnd", sels, value = T)
      ones.reduced <- gsub("Max", "", ones)
      twos.reduced <- gsub("Mnd", "", twos)
      ones <- ones[ones.reduced %in% twos.reduced]
      twos <- twos[twos.reduced %in% ones.reduced]
      invert <- ifelse(input$invert_bone == "inverted", T, F)
      one <- ifelse(invert, "maxilla", "mandible")
      two <- ifelse(invert, "mandible", "maxilla")}
    ages <- unique(substr(sels, 1, 4))
    if (length(ages)) age.1 <- min(ages) else age.1 <- 10.5
    young.col <- paste0("age", age.1)
    if (inputCONFIG$comp == "age"){
      age.2 <- max(ages)
      target.col <- paste0("age", age.2)
      ones <- grep(age.1, sels, value = T)
      twos <- grep(age.2, sels, value = T)
      ones.reduced <- gsub(age.1, "", ones)
      twos.reduced <- gsub(age.2, "", twos)
      ones <- ones[ones.reduced %in% twos.reduced]
      twos <- twos[twos.reduced %in% ones.reduced]
      invert <- ifelse(input$invert_age == "inverted", T, F)
      one <- ifelse(invert, paste0("E", age.2), paste0("E", age.1))
      two <- ifelse(invert, paste0("E", age.1), paste0("E", age.2))}
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
    if (!nrow(design) | !target.col %in% colnames(design) | !all(sels %in% c(ones, twos))) return({
      output$ma.plot <- renderPlot(plot.null())
      output$heatmap <- renderPlot(plot.null())
      data.frame(NULL)})
    output$numsel <- reactive(length(c(ones, twos)))  
    outputOptions(output, "numsel", suspendWhenHidden = FALSE)

#### 
#### HERE is when the data is processed from full raw set to raw selected set
####
    dat.sel <- dat[, c(ones, twos, "symbol")]
    if (input$summary == "AVE"){
      dat.sel <- ddply(dat.sel, .(symbol), numcolwise(mean))
      dat.sel <- dat.sel[dat.sel$symbol != "NA", ]
      rownames(dat.sel) <- dat$probeset[match(dat.sel$symbol, dat$symbol)]}
    if (input$summary == "A"){
      dat.sel$A <- rowMeans(dat.sel)
      dat.sel$P <- rownames(dat.sel)
      dat.single <- dat.sel[dat.sel$symbol %in% single.genes, ]
      dat.multiple <- dat.sel[!dat.sel$symbol %in% single.genes, ]
      dat.agg <- aggregate(A ~ symbol, dat.multiple, max)
      dat.multiple <- merge(dat.agg, dat.multiple, sort = F)
      rownames(dat.multiple) <- dat.multiple$P
      dat.sel <- rbind(dat.single, dat.multiple)
      dat.sel <- dat.sel[, colnames(dat.sel)[!colnames(dat.sel) %in% c("A", "P")]]}
    if (input$summary == "M"){
      dat.sel$MM <- abs(rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones]))
      dat.sel$P <- rownames(dat.sel)
      dat.single <- dat.sel[dat.sel$symbol %in% single.genes, ]
      dat.multiple <- dat.sel[!dat.sel$symbol %in% single.genes, ]
      dat.agg <- aggregate(MM ~ symbol, dat.multiple, max)
      dat.multiple <- merge(dat.agg, dat.multiple, sort = F)
      rownames(dat.multiple) <- dat.multiple$P
      dat.sel <- rbind(dat.single, dat.multiple)
      dat.sel <- dat.sel[, colnames(dat.sel)[!colnames(dat.sel) %in% c("MM", "P")]]}
    fit <- lmFit(dat.sel[, colnames(dat.sel)[colnames(dat.sel) != "symbol"]], design)
    efit <- eBayes(fit)
    dat.sel$A <- rowMeans(dat.sel[, colnames(dat.sel)[colnames(dat.sel) != "symbol"]])
    dat.sel$M <- rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones])
    if (invert) dat.sel$M <- -dat.sel$M
    if (!is.null(dat.sel$symbol))
      dat.sel$symbol[is.na(dat.sel$symbol)] <-  rownames(dat.sel)[is.na(dat.sel$symbol)]
    dat.sel$color <- "black"
## set input$col1 to 'blue' as placeholder
    dat.sel$color[dat.sel$symbol %in% mousecase(input$gene1) | rownames(dat.sel) %in% mousecase(input$gene1)] <- "blue" 

    lfc <- ifelse(input$log, abs(as.numeric(input$fc)), log2(abs(as.numeric(input$fc))))
    top <- topTable(efit, coef = which(colnames(design) == target.col), 
      lfc = lfc, p.value = as.numeric(input$fdr), num = Inf)
    top <- merge(top, dat[, c("probeset", "symbol", "desc")], by.x = 0, by.y = "probeset")
    top.colored <- top[rownames(top) %in% dat.sel$probeset[dat.sel$color != "black"], ]
    if (nrow(top.colored)){
      top.black <- top[1:min(nrow(top), as.numeric(input$max)), ]
      top <- merge(top.colored, top.black, all = T)}
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
      colnames(top) <- c("Probeset", "Gene", "Ave.Expr", "Log2FC", "Fold.Change", "Adj.P.Val", "Description")}
    else top <- data.frame(NULL)
    dat.top <- dat.sel[rownames(dat.sel) %in% top$Probeset, ]

#### 
#### HERE is when the data is now already semi-processed 
####
    output$ma.plotly <- renderPlotly({
      mList <-generateMAplotjson_f(ones,twos,inputCONFIG,dat.sel,dat.top) 
      my.ylab <- paste0(one,
            paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""),
            ifelse(input$log, "Log2 Fold Change", "Fold Change"),
            paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""), two)
      my.xlab <- "Average Expression"
      my.title <- paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " ") 
      generatePlotlyMAplot_f(mList,my.xlab,my.ylab,my.title)
    })

    output$ma.plot <- renderPlot({
      ma.plot <- function(){
        par(mar = c(5, 5, 5, 8))
        lim <- max(abs(range(dat.sel$M)))

        plot(NULL, cex = 0.9,
          xaxt = "n", yaxt = "n",
          xlab = "Average Expression",
          xlim = range(dat.sel$A),
          ylim = c(-lim * 1.05, lim * 1.05), 
          ylab = paste0(one,
            paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""),
            ifelse(input$log, "Log2 Fold Change", "Fold Change"),
            paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""), two),
          cex.lab = 1.2, font.main = 1, cex.main = 1.2,
          main = paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " "))
        points(dat.sel$A[dat.sel$color == "black"], dat.sel$M[dat.sel$color == "black"], pch = 21, cex = 0.9)
        points(dat.sel$A[dat.sel$color != "black"], dat.sel$M[dat.sel$color != "black"], pch = 21, cex = 0.9,
          col = dat.sel$color[dat.sel$color != "black"])
        if (nrow(dat.top)){
          points(dat.top$A, dat.top$M, pch = 21, cex = 0.9, col = dat.top$color, bg = dat.top$color)
          text(dat.top$A, dat.top$M, dat.top$symbol, col = dat.top$color, 
            pos = ifelse(dat.top$M > 0, 3, 1), offset = 0.4, cex = 0.9)}
        ats.0 <- seq(1, 9, 1)
        ats <- c(-1 * rev(ats.0), 0, ats.0)
        log.labs <- ats
        log.labs[log.labs > 0] <- paste0("+", log.labs[log.labs > 0])
        raw.labs <- 2 ^ ats
        raw.labs[raw.labs < 1] <- -1 / raw.labs[raw.labs < 1]
        raw.labs[raw.labs > 1] <- paste0("+", raw.labs[raw.labs > 1])
        axis(1, seq(0, 20, 1))
        if (input$log) axis(2, ats, log.labs)
        else axis(2, ats, raw.labs, cex.axis = ifelse(lim > 6, 0.75, 1))
        abline(h = lfc * c(-1, 1), lty = 2)}
      ma.plot()
    })
      
    output$heatmap <- renderPlot({
      heatmap <- function(){
      withProgress({
        if (nrow(dat.top) < 2) return(plot.null())
        dat.heat <- as.matrix(dat.top[, sels])
        dat.heat <- dat.heat - mean(dat.heat)
        if (input$heatadjust){
          for (col in colnames(control)[colnames(control) != "(Intercept)"]){
            for (v in unique(control[, col])){
              mm <- mean(dat.heat[, rownames(control)[control[, col] == v]])
              dat.heat[, rownames(control)[control[, col] == v]] <- 
                dat.heat[, rownames(control)[control[, col] == v]] - mm}}}
        if (input$heatscale %in% c("MC", "Z")) dat.heat <- sweep(dat.heat, 1, rowMeans(dat.heat))
        if (input$heatscale == "Z") dat.heat <- sweep(dat.heat, 1, apply(dat.heat, 1, sd), "/")
        weights <- (ncol(dat.heat):1) + 10000
        weights[colnames(dat.heat) %in% twos] <- 
          weights[colnames(dat.heat) %in% twos] + ifelse(xor(invert, inputCONFIG$comp == "bone"), -100, 100)
        extreme <- ceiling(10 *  max(dat.heat)) / 10

        heatmap.22(t(dat.heat),
          Rowv = as.integer(weights),
          col = ifelse(input$heatcol == "rg", redgreen, greenred),
          cexCol = min(1.5, 0.2 + 1/log10(nrow(dat.heat))),
          scale = "none", key = T, key.title = NA, lwid = c(1, 4), lhei = c(1.5, 4),
          density.info = "density", densadj = 0.5, denscol = "white",
          key.ylab = NA, key.xlab = ifelse(input$log, "Log2 Fold Change", "Fold Change"),
          key.par = list(cex.lab = 1.25, cex.axis = 1, tcl = -0.35, mgp = c(2, 1, 0)),
          key.ytickfun = function() list(labels = F, tick = F),
          key.xtickfun = function(){
            breaks <- parent.frame()$breaks
            list(at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
               labels = round(seq(-extreme, extreme, length = 6), 1), 
               padj = -0.25)},
          trace = "none", margins = c(5, 10),

          distfun.row = dist,
          distfun.col = function(...) {
            if (input$distfun == "E") return(dist(...))
            if (input$distfun == "AC") return(cor.dist(..., abs = T))
            if (input$distfun == "C") return(cor.dist(..., abs = F))},
          labCol = dat.top$symbol, colCol = dat.top$color)
      incProgress(1)}, value = NULL, message = "Calculations in progress...")}
      heatmap()
    })
      
    output$heatmapPlotly <- renderPlotly({
        if (nrow(dat.top) < 2) return(plot.null())
        dat.heat <- as.matrix(dat.top[, sels])
        dat.heat <- dat.heat - mean(dat.heat)
        if (input$heatadjust){
          for (col in colnames(control)[colnames(control) != "(Intercept)"]){
            for (v in unique(control[, col])){
              mm <- mean(dat.heat[, rownames(control)[control[, col] == v]])
              dat.heat[, rownames(control)[control[, col] == v]] <-
                dat.heat[, rownames(control)[control[, col] == v]] - mm}}}
        if (input$heatscale %in% c("MC", "Z")) dat.heat <- sweep(dat.heat, 1, rowMeans(dat.heat))
        if (input$heatscale == "Z") dat.heat <- sweep(dat.heat, 1, apply(dat.heat, 1, sd), "/")
        weights <- (ncol(dat.heat):1) + 10000
        weights[colnames(dat.heat) %in% twos] <-
          weights[colnames(dat.heat) %in% twos] + ifelse(xor(invert, inputCONFIG$comp == "bone"), -100, 100)
        extreme <- ceiling(10 *  max(dat.heat)) / 10

## row are genes and col are samples
      distfun.col <- dist
      distfun.row <- function(...) {
        if (input$distfun == "AC") return(cor.dist(..., abs = T))
        if (input$distfun == "C") return(cor.dist(..., abs = F))
        return(dist(...))
      }
      my.color <- input$heatcol 
      my.xlab <- ifelse(input$log, "Log2 Fold Change", "Fold Change")
      hList <-generateHeatmapjson_f(ones,twos,inputCONFIG,dat.top,dat.heat) 
      generatePlotlyHeatmap_f(hList, distfun.row, distfun.col,my.color,my.xlab)
    })

    rownames(top) <- NULL
    output$download.table <- downloadHandler("Facebase_Microarray_Table.csv",
      function(file) write.csv(as.data.frame(top), file, row.names = F))
    top$Color <- dat.top$color[match(top$Probeset, rownames(dat.top))]
    top <- formatStyle(datatable(top), "Color", target = "row", 
      color = styleEqual(
        c(0, "blue"),
        c("white", "blue")))
  })
    
}
    
shinyApp(ui = ui, server = server)
