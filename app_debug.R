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

### https://datatables.net/manual/tech-notes/7
### $.fn.dataTable.ext.errMode = 'throw';

eontime <-Sys.time()

mycss <- "
#spinner-container {
  position: relative;
}
#loading-spinner {
  position: absolute;
  left: 50%;
  top: 50%;
  z-index: -1;
  margin-top: -33px;  /* half of the spinner's height */
  margin-left: -33px; /* half of the spinner's width */
}
"

ui <- fluidPage( 
useShinyjs(),
tags$style(HTML(mycss)),
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
h5("eonTime"),
verbatimTextOutput("eonTimeText"),
h5("beginTime"),
verbatimTextOutput("beginTimeText"),
h5("LoadTime"),
verbatimTextOutput("loadTimeText"),
    conditionalPanel( condition="$('html').hasClass('shiny-busy') && !output.datloaded",
        helpText("Shiny is initializing ...")),
    conditionalPanel( condition="!$('html').hasClass('shiny-busy') && !output.datloaded",
        helpText("Data is loading ...")),
    conditionalPanel("output.datloaded",
        helpText("Data is loaded."))
  ),
  tags$form(class = "well",
    bootstrapPage(
      div(style="display:inline-block; width:70%", textInput("gene1", "Highlight one or more genes,", NULL)),
      div(style="display:inline-block;", checkboxInput("go", "go", F)),
#      div(style="display:inline-block;", actionButton(inputId="go", label="on")),
      div(style="padding-bottom:3px", "(eg. Myl2 Rgs5 Tnnt2)")
  )),
  tags$form(class = "well",
    bootstrapPage(
      conditionalPanel("output.comp == 'place'",
        div(tags$label("Select direction for comparing place,")),
        radioButtons("invert_place", br(), 
          c("proximal down, distal up" = "normal",
            "distal down, proximal up" = "inverted"), "normal", inline = T)
      ),
      conditionalPanel("output.comp == 'bone'",
        div(tags$label("Select direction for comparing bone,")),
        radioButtons("invert_bone", br(),
          c("mandible down, maxilla up" = "normal",
            "maxilla down, mandible up" = "inverted"), "normal", inline = T)
      ),
      conditionalPanel("output.comp == 'age'",
        div(tags$label("Select direction for comparing age,")),
        radioButtons("invert_age", br(), 
          c("earliest down, latest up" = "normal",
            "latest down, earliest up" = "inverted"), "normal", inline = T)
      ),
      conditionalPanel("output.numsel > 6",
        checkboxInput("heatadjust", "For heatmap, factor out comparisons other than the one selected", value = T)
      )
  )),
  tags$form(class= "well",
    bootstrapPage(
      div(style="display:inline-block; width:70%", textInput("fc", "Fold change cut-off", 2)),
      div(style="display:inline-block", checkboxInput("log", "log2", T)),
      textInput("fdr", "False discovery rate", 0.01, width = "70%"),
      textInput("max", "Max. DE genes/probesets", Inf, width = "70%"),
      radioButtons("summary", "For each gene, show",
        choices = c( "all probesets" = "Z", "average probeset" = "AVE",
                     "most highly expressed probeset" = "A",
                     "most differentially expressed probeset" = "M"),
                                          selected = "Z", inline = T),
      sliderInput(inputId='bkPts',label = "Show base points in percents,",
                      value=20, min=0, max=100)
  )),
  tags$form(class= "well",
    bootstrapPage(
      radioButtons("heatcol", "For heatmap, high values are",
        choices = c("green" = "rg", "red" = "gr"), selected = "rg", inline = T),
      radioButtons("distfun", "For heatmap, cluster genes using",
        choices = c("Euclidean distance" = "E", "absolute correlation" = "AC",
                    "correlation" = "C"), selected = "E", inline = T),
      radioButtons("heatscale", "For heatmap, scale genes by",
        choices = c( "none" = "N", "mean-centering" = "MC", "Z-score" = "Z"),
                                         selected = "Z", inline = T)
  ))),
  mainPanel(width = 8,

  conditionalPanel("$('html').hasClass('shiny-busy')",
      div(style="height:500px; border:2px solid green", id = "spinner-container",
          tags$img(src = "35.gif", id = "loading-spinner"))
  ),

div(style="border:2px solid blue", plotlyOutput("ma.plotly", height = "500px")),
div(style="border:2px solid blue", plotlyOutput("heatmapPlotly", height = "450px")),

    downloadButton("download.table", "Download table"),  br(), br(),
    DT::dataTableOutput("table", width = "90%")
  )
) ## sidebarLayout
)


server <- function(input, output, session){

begintime <- Sys.time()
output$beginTimeText <- renderText({ paste0(begintime) })
output$eonTimeText <- renderText({ paste0(eontime) })

loginfo("XXX server eontime")
loginfo(eontime)
loginfo("XXX server, begintime")
loginfo(begintime)


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
  source("www/heatmap_util.R")
  source("www/maplot_util.R")

## load the data and setup serverCFG
datasetInput <- reactive({

loginfo("XXX datasetInput -starttime")
    query<- parseQueryString(session$clientData$url_search)

    serverCFG <- processConfig_f(query)

    output$urlText <- renderText({ serverCFG$url })
    con <- url(serverCFG$url)

    loadtime <- system.time(load(con))
output$loadTimeText <- renderText({ paste0(loadtime) })
    output$datloaded <- reactive (is.numeric(nrow(dat)))
    outputOptions(output, "datloaded", suspendWhenHidden = FALSE)
     
    output$selText <- renderText({ paste0(serverCFG$sel) })
    output$compText <- renderText({ paste0(serverCFG$comp) })

    output$comp <- reactive(c(serverCFG$comp))
    outputOptions(output, "comp", suspendWhenHidden = FALSE)

    ones <- serverCFG$ones
    twos <- serverCFG$twos
    output$numsel <- reactive(length(c(ones, twos)))
    outputOptions(output, "numsel", suspendWhenHidden = FALSE)

    rownames(dat) <- dat[, "probeset"]
    dlist <- list( dd=dat, ss=serverCFG )
loginfo("  EEE datasetInput")
    dlist
})

###  force the ball to roll
observe({
loginfo("XXX get the real start...")
    dlist <- isolate(datasetInput())
})


createTop <- reactive({
loginfo("XXX createTop")
    dlist <- isolate(datasetInput())
    dat <- dlist$dd
    serverCFG <- dlist$ss

## get preset 
    ones <- serverCFG$ones
    twos <- serverCFG$twos
    sels <- serverCFG$sels
    target.col <- serverCFG$target.col

    plist <- isolate(createDesignControl())
    design <- plist$dd

    if (!nrow(design) | !target.col %in% colnames(design) | !all(sels %in% c(ones, twos))) {
      return( list() )
    }

# these are what dat.top react to
    serverCFG$log <- input$log
    serverCFG$fc <- input$fc
    serverCFG$fdr <- input$fdr
    serverCFG$max <- input$max
    if(input$go)
      serverCFG$gene1 <- isolate(input$gene1)
      else serverCFG$gene1 <- ""

    if(serverCFG$comp == 'place') {
       serverCFG$inverted <- input$invert_place
    }
    if(serverCFG$comp == 'bone') {
       serverCFG$inverted <- input$invert_bone
    }
    if(serverCFG$comp == 'age') {
       serverCFG$inverted <- input$invert_age
    }

    tlist = makeTop_f(dat,serverCFG,design)
loginfo("   EEE createTop")
    return(tlist)
})

createDesignControl <- reactive({
    dlist <- isolate(datasetInput())
    serverCFG <- dlist$ss
    dat <- dlist$dd
    plist <- makeDesignControl_f(dat,serverCFG)
    return(plist)
})

#http://stackoverflow.com/questions/31920286/effectively-debugging-shiny-apps
output$table <- DT::renderDataTable({
loginfo("XXX table")
    tlist <- createTop()
    top <- tlist$tt
    dat.top <- tlist$dd
### 
    rownames(top) <- NULL
    output$download.table <- downloadHandler("Facebase_Microarray_Table.csv",
      function(file) write.csv(as.data.frame(top), file, row.names = F))
    top$Color <- dat.top$color[match(top$Probeset, rownames(dat.top))]
    top <- formatStyle(datatable(top), "Color", target = "row", 
      color = styleEqual(
        c(0, "blue"),
        c("white", "blue")))
})


#### 
#### HERE is when the data is now already semi-processed 
####
output$ma.plotly <- renderPlotly({
loginfo("XXX ma.plotly")
    dlist <- isolate(datasetInput())
    dat <- dlist$dd
    serverCFG <- dlist$ss

## get preset 
    ones <- serverCFG$ones
    twos <- serverCFG$twos
    one <- serverCFG$one
    two <- serverCFG$two

    tlist <- createTop()
    dat.top <- tlist$dd
    dat.sel <- tlist$ss

    lfc <- ifelse(input$log, abs(as.numeric(input$fc)), log2(abs(as.numeric(input$fc))))
### 
    mList <-generateMAplotjson_f(dat.sel,dat.top) 
    my.ylab <- paste0(one,
        paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""),
        ifelse(input$log, "Log2 Fold Change", "Fold Change"),
        paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""), two)
    my.xlab <- "Expression"
    my.title <- paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " ") 
    my.lfc <- lfc
    my.bkPts <- input$bkPts
    generatePlotlyMAplot_f(mList,my.xlab,my.ylab,my.title,my.lfc,my.bkPts)
})


output$heatmapPlotly <- renderPlotly({
loginfo("XXX heatmap")
    dlist <- isolate(datasetInput())
    serverCFG <- dlist$ss
    tlist <- createTop()
    dat.top <- tlist$dd
    plist <- isolate(createDesignControl())
    control <- plist$cc
    
###
    if (nrow(dat.top) < 2) return(plot.null())

    serverCFG$heatadjust <- input$heatadjust
    serverCFG$heatscale <- input$heatscale

    hlist <- makeHeat_f(dat.top, serverCFG, control)
    dat.heat <- hlist$hh

## row are genes and col are samples
    distfun.col <- dist
    distfun.row <- function(...) {
        if (input$distfun == "AC") return(cor.dist(..., abs = T))
        if (input$distfun == "C") return(cor.dist(..., abs = F))
        return(dist(...))
    }
    my.color <- input$heatcol 
    my.xlab <- ifelse(input$log, "Log2 Fold Change", "Fold Change")
    my.symbols <- dat.top$symbol
    hList <-generateHeatmapjson_f(my.symbols, dat.heat) 
    generatePlotlyHeatmap_f(hList, distfun.row, distfun.col,my.color,my.xlab)
})

}

shinyApp(ui = ui, server = server)

