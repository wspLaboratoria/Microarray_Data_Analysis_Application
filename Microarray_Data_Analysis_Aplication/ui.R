library(shiny)
library(shinyFiles)
library(htmltools)
library(vroom)
library(Biobase)
library(shinycustomloader)


ui <- fluidPage(
    
    titlePanel("Microarray Analysis Aplication"),
    
    sidebarLayout(
        
        sidebarPanel(width = 2,
            
            ## upload text file with phenoData
            fileInput(inputId = "phenoData", 
                      multiple = TRUE, 
                      accept = ".txt",
                      label = "Upload phenoData (.txt) file"),
            
            ## upload .CEL files with microarray data
            fileInput(inputId = "celInfiles", 
                      multiple = TRUE, 
                      accept = ".CEL",
                      label = "Upload the .CEL file list"),
            
            ## selectBox for choosing normalization method
            selectInput("normalizeMethod",
                        label = "Choose normalization method",
                        choices = list("RMA","GCRMA","MAS","VSN"),
                        selected = 1),
            
            ## selectBox for choosing kind of analysis
            selectInput("analysisChoice",
                        label = "Choose type of analysis",
                        choices = list("Differential Gene Expression","C-Means Clustering","Hierarchical Clustering"),
                        selected = 1),
            
        ),
        
        mainPanel(width = 10,
            tabsetPanel(
                tabPanel(title = "Quality Control",
                         icon = icon("check-circle"),
                         width = 1000,
                         
                         
                             fluidRow(
                                 column(6, align = "center", strong("RAW DATA")),
                                 column(6, align = "center", strong("NORMALIZED DATA"))
                             ),
                                 
                             fluidRow(
                                 column(6,withLoader(plotOutput("boxplotRaw")),type = "html",loader="dmaspin"),
                                 column(6,withLoader(plotOutput("boxplotNorm")),type = "html",loader="dmaspin")
                             ),
                             fluidRow(
                                 column(6,withLoader(plotOutput("histRaw")),type = "html",loader="dmaspin"),
                                 column(6,withLoader(plotOutput("histNorm")),type = "html",loader="dmaspin")
                             ),
                             fluidRow(
                                 column(6,withLoader(plotOutput("maPlotRaw")),type = "html",loader="dmaspin"),
                                 column(6,withLoader(plotOutput("maPlotNorm")),type = "html",loader="dmaspin")
                            )       
                        
                        ),
                
                tabPanel(title = uiOutput("tabTitle"),
                         plotOutput(""))
                
            )
        )
    )
)