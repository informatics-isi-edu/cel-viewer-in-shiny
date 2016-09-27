eontime <-Sys.time()
library(bioDist)
library(limma)
library(DT)
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

eon2time <-Sys.time()

ui <- fluidPage( 
useShinyjs(),
titlePanel("Gene Expression Comparison Shiny App"),
sidebarLayout(
  div(class = "col-sm-4",
  tags$form(class = "well",
    h5("URL"),
    verbatimTextOutput("urlText"),
    h5("Selection"),
    verbatimTextOutput("selText"),
    h5("Comparison"),
    verbatimTextOutput("compText"),
h5("LoadTime"),
verbatimTextOutput("loadTimeText"),
h5("eonTime"),
verbatimTextOutput("eonTimeText"),
h5("eon2Time"),
verbatimTextOutput("eon2TimeText"),
h5("beginTime"),
verbatimTextOutput("beginTimeText"),
h5("startTime"),
verbatimTextOutput("startTimeText"),
h5("preloadTime"),
verbatimTextOutput("preloadTimeText"),
h5("endTime"),
verbatimTextOutput("endTimeText"),
    conditionalPanel("!output.datloaded",
        helpText("Data is loading (may take a minute) ...")),
    conditionalPanel("output.datloaded", 
        helpText("Data is loaded."))),

   tags$form(class = "well",
       bootstrapPage(
        div(tags$label("Highlight one or more genes, ")),
        div(style="display:inline-block; width:70%", textInput("gene1", NULL, NULL)),
        div(style="padding-bottom:3px", "(eg. Myl2 Rgs5 Tnnt2)")
       )),
  tags$form(class = "well",
    conditionalPanel("output.comp == 'place'",
     div(tags$label("Select direction for comparing place,"))),
    conditionalPanel("output.comp == 'bone'",
     div(tags$label("Select direction for comparing bone,"))),
    conditionalPanel("output.comp == 'age'",
     div(tags$label("Select direction for comparing age,"))),
  bootstrapPage(
    conditionalPanel("output.comp == 'place'",
     radioButtons("invert_place", br(), 
         c("proximal down, distal up" = "normal",
           "distal down, proximal up" = "inverted"), "normal", inline = T)
    ),
    conditionalPanel("output.comp == 'bone'",
       radioButtons("invert_bone", br(),
           c("mandible down, maxilla up" = "normal",
             "maxilla down, mandible up" = "inverted"), "normal", inline = T)
    ),
    conditionalPanel("output.comp == 'age'",
       radioButtons("invert_age", br(), 
           c("earliest down, latest up" = "normal",
             "latest down, earliest up" = "inverted"), "normal", inline = T)
    )
   ),
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

    div(style="border:2px solid blue", plotlyOutput("ma.plotly", height = "500px")),

    div(style="border:2px solid blue", plotlyOutput("heatmapPlotly", height = "450px")),

    downloadButton("download.table", "Download table"),  br(), br(),
    DT::dataTableOutput("table", width = "90%"))
) ## sidebarLayout
)


server <- function(input, output, session){
begintime <- Sys.time()
output$beginTimeText <- renderText({ paste0(begintime) })


#http://deanattali.com/blog/advanced-shiny-tips/#shinyjs

session$onSessionEnded(stopApp)

printLogJs <- function(x, ...) { 
                               logjs(x) 
                               T
                               }
addHandler(printLogJs)

pdf(NULL)

## process args that were passed in via url

##https://dev.facebase.org/shiny/apps/url/?url=hello-world
##
##https://dev.facebase.org/shiny/apps/usc/?url=https://dev.facebase.org/~mei/data/urlexprs.R
  source("www/util.R")
  source("www/old_util.R")
  source("www/heatmap_util.R")
  source("www/maplot_util.R")

  age <- factor(rep(c(10.5, 11.5, 12.5, 13.5, 14.5), each = 12))
  bone <- factor(rep(rep(c("Max", "Mnd"), each = 6), 5))
  place <- relevel(factor(rep(rep(c("D", "P"), each = 3), 10)), "P")
  full.design <- model.matrix(~age + bone + place)
  rownames(full.design) <- paste0(age, bone, place, 1:3)

serverCFG <- list( sel=c("E10.5_Mnd_D","E10.5_Mnd_P"), comp= c("place"), url= NULL)

  processForData = reactive( {

## process input config part
starttime <- Sys.time()
output$startTimeText <- renderText({ paste0(starttime) })

  query<- parseQueryString(session$clientData$url_search)

  if (!is.null(query[['url']])) {
    loginfo("AAA 1")
    serverCFG$url <- query[['url']]
    output$urlText <- renderText({ paste0(serverCFG$url) })
    con <- url(serverCFG$url)

    preloadtime <- Sys.time()
output$preloadTimeText <- renderText({ paste0(preloadtime) })
    loadtime <- system.time(load(con))
output$loadTimeText <- renderText({ paste0(loadtime) })
    loginfo("AAA 2")
  } 
  loginfo("AAA 3")

  if (!is.null(query[['sel']]))  {
# start with "E10.5_Mnd_D E10.5_Mnd_P"
# into  "E10.5_Mnd_D" "E10.5_Mnd_P"
    s <- query[['sel']]
    s <- gsub("\"", "", s)
    ss <- strsplit(s," ")
    ss <- ss[[1]]
    t <- character()
    for(s in ss) t <- c(t, s)
    serverCFG$sel <- t
    output$selText <- renderText({ paste0(serverCFG$sel) })
  }
  loginfo("AAA 4")

  if (!is.null(query[['comp']])) { 
    s <- query[['comp']]
    s <- gsub("\"", "", s)
    serverCFG$comp <- s
    output$compText <- renderText({ paste0(serverCFG$comp) })
  }
  loginfo("AAA 5")
  dat
  })
  loginfo("AAA 6")
  dat <- isolate(processForData())
  loginfo("AAA 7")

  rownames(dat) <- dat[, "probeset"]
  genes.tab <- xtabs(~unlist(dat$symbol))
  single.genes <- names(genes.tab)[genes.tab == 1]

  output$datloaded <- reactive(is.numeric(nrow(dat)))
  outputOptions(output, "datloaded", suspendWhenHidden = FALSE)
  output$comp <- reactive(serverCFG$comp)
  outputOptions(output, "comp", suspendWhenHidden = FALSE)

  loginfo("AAA 8")
  
#http://stackoverflow.com/questions/31920286/effectively-debugging-shiny-apps
  output$table <- DT::renderDataTable({
    loginfo("XXX renderDataTAble")
    sel <- serverCFG$sel
    sel <- substr(sel, 2, nchar(sel))
    sel <- gsub("_", "", sel)
    sels <- character()
    for (s in sel) sels <- c(sels, paste0(s, 1:3))
    if (serverCFG$comp == "place"){ 
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
    if (serverCFG$comp == "bone"){
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
    if (serverCFG$comp == "age"){
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
      output$ma.plotly <- renderPlotly(plot.null())
      output$heatmapPloty <- renderPlotly(plot.null())
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
## initialize that the gene1 to be empty
    gene1 <- c(strsplit(input$gene1," "))
    gene1 <- gene1[[1]]
    if( !(is.null(gene1)) ) {
      for (g in gene1) {
        dat.sel$color[dat.sel$symbol == mousecase(g) | rownames(dat.sel) == mousecase(g)] <- "blue"
      }
    }

    lfc <- ifelse(input$log, abs(as.numeric(input$fc)), log2(abs(as.numeric(input$fc))))

    top <- topTable(efit, coef = which(colnames(design) == target.col), 
      lfc = lfc, p.value = as.numeric(input$fdr), num = Inf)
    top <- merge(top, dat[, c("probeset", "symbol", "desc")], by.x = 0, by.y = "probeset")
    top.colored <- top[rownames(top) %in% dat.sel$probeset[dat.sel$color != "black"], ]
## only take the max into account if there is a highlighted/colored top entry
    if (nrow(top.colored)){
loginfo("BBB")
      top.black <- top[1:min(nrow(top), as.numeric(input$max)), ]
      top <- merge(top.colored, top.black, all = T)}
loginfo("CCC")
loginfo(nrow(top))
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
        loginfo("XXX ma.plotly")
        mList <-generateMAplotjson_f(ones,twos,serverCFG,dat.sel,dat.top) 
        my.ylab <- paste0(one,
            paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""),
            ifelse(input$log, "Log2 Fold Change", "Fold Change"),
            paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""), two)
        my.xlab <- "Average Expression"
        my.title <- paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " ") 
        my.lfc <- lfc
      generatePlotlyMAplot_f(mList,my.xlab,my.ylab,my.title,my.lfc)
    })

    output$heatmapPlotly <- renderPlotly({
        loginfo("XXX heatmap")
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
          weights[colnames(dat.heat) %in% twos] + ifelse(xor(invert, serverCFG$comp == "bone"), -100, 100)
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
      hList <-generateHeatmapjson_f(ones,twos,serverCFG,dat.top,dat.heat) 
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

endtime <- Sys.time()
output$endTimeText <- renderText({ paste0(endtime) })

output$eonTimeText <- renderText({ paste0(eontime) })
output$eon2TimeText <- renderText({ paste0(eon2time) })
#  }) ## outermost observe
    
}
    
shinyApp(ui = ui, server = server)
