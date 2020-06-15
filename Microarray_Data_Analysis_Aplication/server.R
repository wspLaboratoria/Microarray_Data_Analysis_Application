# shiny packages
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(htmltools)

# normalization of data
library(affy)
library(gcrma) 
library(limma)
library(oligo)

# to create plotly plots
library(plotly)
library(ggplot2)

# ### hierarhical clastering - hc.R ###
library(dendextend)
library(heatmaply)
library(factoextra)
library(ggdendro)
library(colorspace)

### C-means clastering - kmen.R ###
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra)
library(amap)

server <- function(input, output, session) {
    
    ### setting paths and options ###
    options(shiny.maxRequestSize=30*1024^2) 
    
    ###  SOURCE FOR FUNCTIONS ###
    ## Quality Control 
    source(paste(getwd(),"/R/readData.R", sep=""))
    source(paste(getwd(),"/R/normalization.R", sep=""))
    source(paste(getwd(),"/R/boxplot.R", sep=""))
    source(paste(getwd(),"/R/histPlot.R", sep=""))
    source(paste(getwd(),"/R/maPlot.R", sep=""))
    
    
    ## Differential Gene Expression ###
    source(paste(getwd(),"/R/pca_genes.R", sep=""))
    source(paste(getwd(),"/R/hitmap_genes.R", sep=""))
    source(paste(getwd(),"/R/dataTableDifferencialGenes.R", sep=""))
    source(paste(getwd(),"/R/dataTablesignalTransductionPathway.R", sep=""))
    
    ### Hierarchical Clustering ###
    source(paste(getwd(),"/R/hc.R", sep=""))
    
    ### C-Means Clustering ###
    source(paste(getwd(),"/R/kmen.R", sep=""))
   
    ### reactive variables ###
    dataUploaded <- reactive({
        validate(
            need(input$celInfiles, 'Please enter microarray data (.CEL files)'),
            need(input$phenoData, 'Please enter pheno data (.txt file)')
        )
        (!(is.null(input$phenoData) && is.null(input$celInfiles)))
    })
    
    dataUploadedNoComment <- reactive({
        validate(
            need(input$celInfiles, ''),
            need(input$phenoData, '')
        )
        return(!(is.null(input$phenoData) && is.null(input$celInfiles)))
    })
    
    resultKmenCluster <- reactive({
        if (dataUploadedNoComment()) {
            kmen(dataNorm(), input$kmen_numGenes, input$kmen_numCluster, input$kmen_distMethod)
        }
    })
    
    resultHCluster <- reactive({
        if (dataUploadedNoComment()) {
            hc(dataNorm(), input$h_numGenes, input$h_distMethod, input$h_clusterMethod, input$h_numCluster)
        }
    })
    
    sampleInData <- reactive({
        data()@phenoData@data$Sample
    })
    
    classInData <- reactive({
        dataNorm()@phenoData@data$CLASS
    })
    
    ### ui parts ###
    output$tabTitle <- renderText({input$analysisChoice})
    
    output$headers <- renderUI({
      if (dataUploadedNoComment()) {
        dashboardBody(
          fluidRow(
            column(6, align = "center", strong("Raw Data")),
            column(6, align = "center", strong("Normalized Data"))
          )
        )
      }
    })

    ### Quality Control ###
    # read data #
    data <- reactive({
        if (dataUploadedNoComment()) {
            readData(input$celInfiles$datapath,input$phenoData$datapath)
        }
    })

    # normalization data #
    dataNorm <- reactive({
        if (dataUploadedNoComment()) {
            normalization(data(),input$normalizeMethod)
        }
    })

    # number of genes for Quality Control #
    output$qc_numGenes <- renderUI({
      if (dataUploadedNoComment()) {
          dashboardBody(
            sliderInput("qc_numGenes", "Number of genes", min = 2, max = dim(dataNorm())[1], value = 100, step = 1, width = '75%')
          )
        }
    })
    outputOptions(output, "qc_numGenes", priority = 10)

    # boxplot #
    output$boxplotRaw <- renderPlotly({
        # use validate with comment to inform user to load files
        if (dataUploaded()) {
          boxplot(data(),dataNorm(),input$qc_numGenes)
        }
    })


    
    output$boxplotNorm <- renderPlotly({
        if (dataUploadedNoComment()) {
            boxplot(dataNorm(),dataNorm(),input$qc_numGenes)
        }
    })

    # histogram #
    output$histRaw <- renderPlotly({
        if (dataUploadedNoComment()) {
            histPlot(data(),dataNorm(),input$qc_numGenes)
        }
    })

    output$histNorm <- renderPlotly({
        if (dataUploadedNoComment()) {
            histPlot(dataNorm(),dataNorm(),input$qc_numGenes)
        }
    })
    
    # maPlot - choosing the sample for plotting #
    output$sampleNamesSelectBox <- renderUI({
        if (dataUploadedNoComment()) {
            
            dashboardBody(
                selectInput("sampleMaPlot",
                            label = "Choose sample for which You want to display maPlot",
                            choices =  sampleInData(),
                            selected = 1,
                            selectize = FALSE)
            )
        }
    })

    output$prevBin <- renderUI({
        if (dataUploadedNoComment()) {
            actionButton("prevBin", label = "Previous")
        }
    })
    
    output$nextBin <- renderUI({
        if (dataUploadedNoComment()) {
            actionButton("nextBin", 
                         label = "Next")
        }
    })
    
    observeEvent(input$prevBin, {
        current <- which(sampleInData() == input$sampleMaPlot)
        if(current > 1){
            updateSelectInput(session, "sampleMaPlot",
                              choices = sampleInData(),
                              selected = sampleInData()[current - 1])
        }
    })
    observeEvent(input$nextBin, {
        current <- which(sampleInData() == input$sampleMaPlot)
        if(current < length(sampleInData())){
            updateSelectInput(session, "sampleMaPlot",
                              choices = sampleInData(),
                              selected = sampleInData()[current + 1])
        }
    })
    
    # maPlot #
    output$maPlotRaw <- renderPlot({
        if (dataUploadedNoComment()) {
            maPlot(data(),input$sampleMaPlot)
        }
    })

    output$maPlotNorm <- renderPlot({
        if (dataUploadedNoComment()) {
            maPlot(dataNorm(),input$sampleMaPlot)
        }
    })

    ### Differential Gene Expression ###
    # TODO: Ola's hitmap
    output$hitmap_genes <- renderPlot({
        # if (dataUploadedNoComment()) {
            # hitmap_genes(dataNorm())
        # }
        hitmap_genes(2)
    })  
    
    output$pca_genes <- renderPlot({
        if (dataUploadedNoComment()) {
            pca_genes(dataNorm())
        }
    })  
    
    output$dataTableDifferencialGenes <- renderTable({
        if (dataUploadedNoComment()) {
            dataTableDifferencialGenes(dataNorm(),input$group1_difGenes,input$group2_difGenes)
        }
    })

    # TODO: Milena's function!
    output$dataTableSignalTransductionPathway <- renderTable({
        # if (dataUploadedNoComment()) {
        # dataTableSignalTransductionPathway(dataNorm())
        # }
        dataTablesignalTransductionPathway(2)
    })  
    
    ### Hierarchical Clustering ###
    output$h_cluster <- renderPlot({
        if (dataUploadedNoComment()) {
            resultHCluster()[[1]]
        }
    })  

    output$hitmap_dendogram_h_cluster <- renderPlotly({
        if (dataUploadedNoComment()) {
            resultHCluster()[[2]]
        }
    })  
    
    ### C-Means Clustering ###
    output$kmen_cluster <- renderPlot({
        if (dataUploadedNoComment()) {
            resultKmenCluster()[[1]]
        }
    })  
    
    
    output$kmen_table <- renderTable({
        if (dataUploadedNoComment()) {
            resultKmenCluster()[[2]]
        }
    })  
    
    ### Layout of the second tabPanel ###
    output$layoutTabPanel <- renderUI({
        typeOfAnalysis <- input$analysisChoice
        if (typeOfAnalysis == "Differential Gene Expression") {
           
            if (dataUploaded()) {
                dashboardBody(
                    
                    fluidRow(
                        column(12, align = "left", strong("Choose two groups to perform t test:"))
                    ),
                    
                    fluidRow(
                        column(3, align = "left", selectInput("group1_difGenes",label = "Group 1",
                                                                choices = unique(classInData()),
                                                                selected = unique(classInData())[[1]],
                                                                selectize = FALSE)),
                        column(3, align = "left", selectInput("group2_difGenes",label = "Group 2",
                                                                choices = unique(classInData()),
                                                                selected = unique(classInData())[[2]],
                                                                selectize = FALSE))
                    ),
                    
                    fluidRow(
                        column(6, align = "center", plotOutput("pca_genes")),
                        
                        # TODO: waiting for real function for Heatmap
                        # column(4,withLoader(plotOutput("hitmap_genes")),type = "html",loader="dmaspin")
                        column(6, align = "center", plotOutput("hitmap_genes"))
                    ),
                    
                    fluidRow(
                        column(6, align = "center", strong("List of differencial genes")),
                        column(6, align = "center", strong("Signal transduction pathway analysis results"))
                    ),
                            
                    fluidRow(
                        column(6, align = "center", tableOutput("dataTableDifferencialGenes")),
                        column(6, align = "center", tableOutput("dataTablesignalTransductionPathway"))
                    )
      
                )
            }
            
        } else if (typeOfAnalysis == "C-Means Clustering") {
            
            if (dataUploaded()) {
                dashboardBody(
                    
                    fluidRow(
                        column(3, align = "center", sliderInput("kmen_numCluster", "Number of clusters:", min = 0, max = dim(dataNorm())[2], 
                                                                value = (dim(dataNorm())[2])/2, step=1)),
                        column(3, align = "center", selectInput("kmen_distMethod",label = "Choose type of distance method",
                                                                choices = list("euclidean", "maximum", "manhattan", "canberra"),
                                                                selected = 1, selectize = FALSE)),
                        column(6, align = "center", sliderInput("kmen_numGenes", "Number of genes:", min = 0, max = dim(dataNorm())[1], value = 100, 
                                                                step = 1, width = '100%'))
                    ),
                    
                  
                    fluidRow(
                        column(6, align = "center", plotOutput("kmen_cluster")),
                        column(6, align = "center", tableOutput("kmen_table"))
                    )
                )
            }
            
        # Hierarchical Clustering
        } else { 
            
            if (dataUploaded()){
                dashboardBody(
    
                    fluidRow(
                        column(4, align = "center", sliderInput("h_numCluster", "Number of clusters:", min = 0, max = dim(dataNorm())[2], 
                                                                value = (dim(dataNorm())[2])/2, step=1)),
                        
                        column(4, align = "center", selectInput("h_distMethod",label = "Choose type of distance method",
                                                                choices = list("euclidean", "maximum", "manhattan", "canberra"),
                                                                selected = 1,selectize = FALSE)),
                        column(4, align = "center", selectInput("h_clusterMethod",label = "Choose type of clustering method",
                                                                choices = list("average","complete","single","centroid"),
                                                                selected = 1, selectize = FALSE)),
                    ),
                    
                    fluidRow(
                        column(12, align = "center", sliderInput("h_numGenes", "Number of genes:", min = 0, max = dim(dataNorm())[1], 
                                                                 value = 100, step = 1, width = '50%'))
                    ),
                    
    
                    fluidRow(
                        column(6, align = "center", plotOutput("h_cluster")),
                        column(6, align = "center", plotlyOutput("hitmap_dendogram_h_cluster"))
                    )
                )
            }
        }
        
    })
}
