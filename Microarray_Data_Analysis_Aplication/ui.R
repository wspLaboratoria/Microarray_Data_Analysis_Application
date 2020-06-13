# shiny package to add loading graphics
library(shinycustomloader)

ui <- fluidPage(
    
    titlePanel("Microarray Analysis Application"),
    
    sidebarLayout(
        
        sidebarPanel(width = 3,
            
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
                        choices = list("RMA","GCRMA","MAS"),
                        selected = 1,
                        selectize = FALSE),
            
            ## selectBox for choosing kind of analysis
            selectInput("analysisChoice",
                        label = "Choose type of analysis",
                        choices = list("Differential Gene Expression","C-Means Clustering","Hierarchical Clustering"),
                        selected = 1,
                        selectize = FALSE),
            
        ),
        
        mainPanel(width = 9,
            tabsetPanel(
                tabPanel(title = "Quality Control",
                         icon = icon("check-circle"),
                         
                         
                             fluidRow(
                                 column(6, align = "center", strong("Raw Data")),
                                 column(6, align = "center", strong("Normalized Data"))
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
                                 column(1, uiOutput("prevBin")),
                                 column(2, uiOutput("sampleNamesSelectBox")),
                                 column(1, uiOutput("nextBin"))
                             ),     
                         
                             fluidRow(
                                 column(6,withLoader(plotOutput("maPlotRaw")),type = "html",loader="dmaspin"),
                                 column(6,withLoader(plotOutput("maPlotNorm")),type = "html",loader="dmaspin")
                            )       
                        
                        ),
                
                
                tabPanel(title = uiOutput("tabTitle"),
                         width = 1000,
                         uiOutput("layoutTabPanel")
                        )
            )
        )
    )
)