# @Yong Kee
# source("global.R")
# source("global2.R")
# library(shinyFiles)
# library(shiny)

#' A user friendly GUI client for oneSENSE package
#'
#' This GUI provides an easy way for Flow Cytometry data analysis using the oneSENSE package. Main parameters for running 'oneSENSE' were integrated in this GUI, and analysis results are launched in Rstudio after submission.
#' Title
#'
#' @return GUI for onesense analysis
#' @export
#' @import shiny
#' @importFrom graphics plot
#' @importFrom utils write.csv
#' @examples
#'
#' if (interactive()) oneSENSE::onesense_GUI()
onesense_GUI <- function() {

ui <- fluidPage(
  titlePanel("OneSENSE"),
  hr(),
  sidebarLayout(
    sidebarPanel(
      ## SELECT FCS RAW DIRECTORY ##
      shinyDirButton('directory', 'Folder select', 'Please select a folder'),
      tags$h4('FCS Folder Selected:'),
      verbatimTextOutput('directorypath'),

      #shinyFilesButton('chooseFCS', 'Choose FCS', 'Please select the FCS files', multiple = TRUE),
      ## SELECT FCS FILES ##

      ## CATERGORY MARKER SELECTION ##
      actionButton("firstMarkers", 'Display Category Markers'),
      actionButton("updatenamescsv", 'Confirm Marker Selection'),

      ## MERGE METHODS ##
      #checkboxGroupInput("merge", "Merge Method", choices = c("All", "Min", "Fixed", "Ceil")),

      ## MERGE INPUT ##
      numericInput("num", "Ceiling Input", value = 500),

      ## BINS ##
      numericInput("bins", "Bins", value = 250),

      ## FREQUENCY HEATPLOT ##
      #checkboxInput("freq", "Do frequency heatplot", value = FALSE),

      ## Submit 1 ##
      actionButton("submit", "Submit"),

      actionButton("submit2", "Coordinate Selection"),
      br(),
      br(),
      actionButton("submit3", "Generate Frequency heatmap")
      ##testoutput##
      #textOutput("text1")
    ),

    mainPanel(tabsetPanel(id = "Tabset1",
      tabPanel("Marker selection", value ="MkrSelect",
               flowLayout(
               checkboxGroupInput("markers1", "Select 1st category markers", ""),
               checkboxGroupInput("markers2", "Select 2nd category markers", ""),
               checkboxGroupInput("markers3", "Select 3rd category markers (Optional)", ""))
      ),
      tabPanel("Coordinate Selection", value = "CoordSelect",
               sidebarLayout(
                 sidebarPanel(
                   uiOutput("markerXdropdown"),
                   uiOutput("markerYdropdown"),
                   tableOutput('table')
                 ),

                 mainPanel(
                   plotOutput("plot", click = "plot_click"),
                   actionButton("writecsv", "Generate CSV")

                 )
               )
               )

      )
    )
      )
)

server <- function(input, output, session) {

  ## SELECT FCS RAW DIRECTORY ##
  # if (.Platform$OS.type == "windows") {
  # volumes <<- c('Root' , "C:/")
  # }
  # else if (.Platform$OS.type == "unix") {
  #   volumes <<- c('Root' , "/root")
  # }
  # else {
  #   volumes <<-c('Root' , "/")
  #   }
  if (.Platform$OS.type == "windows") volumes <- c('Root'="C:/")

  ##need to find out all the platform specific root paths
  shinyDirChoose(input, 'directory', roots=volumes, session=session)
  #shinyFileChoose(input, 'chooseFCS', roots= volumes, session =session)
  output$directorypath <- renderPrint({parseDirPath(roots =volumes, input$directory)})

  ## SELECT FCS FILES ##

  ## SELECT CATEGORY MARKERS ##
  observeEvent(input$firstMarkers, {
    getParameters(parseDirPath(volumes, input$directory))
    updateCheckboxGroupInput(session, "markers1", choices = markers)
    updateCheckboxGroupInput(session, "markers2", choices = markers)
    updateCheckboxGroupInput(session, "markers3", choices = markers)
  })


  ## CONFIRM CATEGORY MARKERS TO WRITE CSV ##
  observeEvent(input$updatenamescsv, {

    paras1 <- c(input$markers1)
    input1 <- rep("", length(markers_desc))
    input1[match(paras1, markers)] <- "Y"
    names1 <- data.frame(Marker=markers_desc, input = input1)

    paras2 <- c(input$markers2)
    input2 <- rep("", length(markers_desc))
    input2[match(paras2, markers)] <- "Y"
    names2 <- data.frame(Marker=markers_desc, input = input2)

    #message(input$markers3)
    paras3 <- c(input$markers3)
    input3 <- rep("", length(markers_desc))
    input3[match(paras3, markers)] <- "Y"
    names3 <- data.frame(Marker=markers_desc, input = input3)

    if(is.null(input$markers3)) {
      final <- cbind(names1, input2 = names2[,2])
      write.csv(final, paste(dirname(parseDirPath(volumes, input$directory)), "names.csv", sep=.Platform$file.sep), row.names = FALSE)
    }
    else {
    final <- cbind(names1, input2 = names2[,2], input3 =names3[,2])
    write.csv(final, paste(dirname(parseDirPath(volumes, input$directory)), "names.csv", sep=.Platform$file.sep), row.names = FALSE)
    }
  })

  ## DO TSNE AND ONESENSE ##
  observeEvent(
    input$submit, {
      FCStSNE(LoaderPATH =parseDirPath(volumes, input$directory),
                           ceil = input$num,
                           FNnames="names.csv",
                           OutputSuffix = "Out",
                           DotSNE = TRUE,
                           DoOneSENSE = TRUE,
                           Bins = input$bins)
      OneSmapperFlour(LoaderPATH =paste(parseDirPath(volumes, input$directory), "_Out", sep =""),
                        Bins = input$bins,
                        doCoords=FALSE, doFreq=FALSE)


      }

  )
  observeEvent(input$submit2, {
    OneSmapperFreq1(LoaderPATH = paste(parseDirPath(volumes, input$directory), "_Out", sep =""))
    getCoords(LoaderPATH = paste(parseDirPath(volumes, input$directory), "_Out", sep =""), FFdata = FFdata)
    updateTabsetPanel(session, "Tabset1", selected = "CoordSelect")

    values <-reactiveValues()
    values$df <- data.frame(Mid, row.names = gckeepnames)

    output$markerXdropdown <- renderUI({
      mydata = gckeepnames
      selectInput("markerX", "Marker X", mydata)
    })
    output$markerYdropdown <- renderUI({
      mydata = gckeepnames
      selectInput("markerY", "Marker Y", mydata)
    })
    observe({ if (!is.null(input$plot_click$x)) {
      isolate(values$df[input$markerX,] <- input$plot_click$x)
    }})

    output$plot<- renderPlot({
      ExpX <- gcdata1[,input$markerX]
      ExpY <- gcdata1[,input$markerY]
      plot(ExpX, ExpY, pch=".", xlab=input$markerX, ylab=input$markerY, cex=0.1, ylim = c(-1,4.5), xlim = c(-1,4.5))
    }, width = 400, height = 400)


    output$table <- renderTable( values$df, rownames = TRUE, digits = 10)

    observeEvent(input$writecsv, {
      write.csv(values$df, paste(dirname(parseDirPath(volumes, input$directory)), "Coords.csv", sep=.Platform$file.sep))
      #stopApp()
    })

  })
  observeEvent(input$submit3, {
    OneSmapperFreq2(LoaderPATH = paste(parseDirPath(volumes, input$directory), "_Out", sep =""), Bins = input$bins)
  })


}



startShinyApp <- shinyApp(ui=ui, server = server)
runApp(startShinyApp)
}
