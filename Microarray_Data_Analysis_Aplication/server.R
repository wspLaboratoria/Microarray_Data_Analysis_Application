library(shiny)
library(shinyFiles)
library(htmltools)
library(vroom)
library(affy);library(gcrma); library(limma)


server <- function(input, output, session) {
    
    ## setting paths and options
    options(shiny.maxRequestSize=30*1024^2) 
    
    ## source for functions 
    source(paste(getwd(),"/R/readData.R", sep=""))
    source(paste(getwd(),"/R/normalization.R", sep=""))
    source(paste(getwd(),"/R/createBoxplot.R", sep=""))
    source(paste(getwd(),"/R/histPlot.R", sep=""))
    source(paste(getwd(),"/R/maPlot.R", sep=""))
    
    output$tabTitle <- renderText({input$analysisChoice})
    
    dataUploaded <- reactive({
        validate(
            need(input$celInfiles, 'Please enter microarray data (.CEL files)'),
            need(input$phenoData, 'Please enter pheno data (.txt file)')
        )
        (!(is.null(input$phenoData) && is.null(input$celInfiles))) 
        
    })
    
    data <- reactive({
        if (dataUploaded()) {
            readData(input$celInfiles$datapath,input$phenoData$datapath)
        }
    })
    
    dataNorm <- reactive({
        if (dataUploaded()) {
            normalization(data(),input$normalizeMethod)
        }
    })
    
    output$boxplotRaw <- renderPlot({
        if (dataUploaded()) {
            createBoxplot(data())
        }
    })
    
    output$boxplotNorm <- renderPlot({
        if (dataUploaded()) {
            createBoxplot(dataNorm())
        }
    })
    
    output$histRaw <- renderPlot({
        if (dataUploaded()) {
            histPlot(data())
        }
    })
    
    output$histNorm <- renderPlot({
        if (dataUploaded()) {
            histPlot(dataNorm())
        }
    })
    
    output$maPlotRaw <- renderPlot({
        if (dataUploaded()) {
            maPlot(data())
        }
    })
    
    output$maPlotNorm <- renderPlot({
        if (dataUploaded()) {
            maPlot(dataNorm())
        }
    })
    
}
