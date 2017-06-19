# Evan Newell 2015
# Edited by Tan Yong Kee


#' Median heatplot generation
#'
#' This returns the median heatplot for each category of markers chosen. If frequency heatplot is selected in the onesense GUI, then both median and frequency heatplots are generated as PDF files in working directory
#'
#' @param LoaderPATH Name of the output file containing fcs files generated from FCStSNE2.R
#' @param Bins Number of bins to sort cells into corresponding heatplot
#' @param doCoords a boolean that allows for frequency heatplot generation
#' @param doFreq a boolean to allow for the frequency heatplot generation. TRUE to run, FALSE to not run.
#'
#' @return PNG files of combined oneSENSE and heatplot.
#'
#'
#' @import webshot
#' @importFrom flowCore exprs read.FCS
#' @importFrom plotly plot_ly subplot export layout
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom stats median
#' @importFrom graphics hist locator
#' @importFrom utils browseURL read.csv write.csv
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom shinyFiles parseDirPath shinyDirButton shinyDirChoose
#'
#'
#' @examples
#' #fcsoutpath <- paste(getwd(), "fcs_Out", sep = .Platform$file.sep)
#' #remove hash symbol to run
#' #OneSmapperFlour(LoaderPATH=fcsoutpath) #remove hash symbol to run
OneSmapperFlour <- function(LoaderPATH ="fcs_Out",
                            Bins = 250,
                            doCoords=FALSE, doFreq=FALSE)

{
  message("Running OneSmapperFlour")
  fs <- read.flowSet(path = LoaderPATH, pattern = '.fcs$') #3 sample Fcs files read
  FcsFileNames <- rownames(keyword(fs, "FILENAME")) #store fcs file names to variable fcsfilenames
  FNumBC <- length(fs) #length of fs = 3
  FFdata <- NULL #Flowframedata
  for (FFs in 1:FNumBC) {
    FFa <- exprs(fs[[FFs]]) #Store data from each fcs file to FFa

    #Fixup column names
    colnames(FFa) <- fs[[FFs]]@parameters$desc #assign column names to each column
    empties <- which(is.na(colnames(FFa)) | colnames(FFa)== " ") #find which columns are NA or space
    colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties] #add the corresponding column names to FFa which are initially empty
    fs[[FFs]]@parameters$desc <- colnames(FFa) #add colnames of FFa to the fcs file
    #fs[[FFs]]@parameters$name <- colnames(FFa)
    #Add file label
    FFa <- cbind(FFa,rep(FFs,dim(FFa)[1])) #rep is replicate, it binds the number associated with the file to the last column
    colnames(FFa)[dim(FFa)[2]] <- "InFile" #give the last column the location of the file
    #Concatenate
    FFdata <- rbind(FFdata,FFa) #adds row data from FFa to FFdata
  }


  lgcl <- logicleTransform(w=0.05, t=16409, m=4.5, a=0)
  ilgcl <- inverseLogicleTransform(trans = lgcl)

  Xx1DtSNEmat <- NULL
  keeptable <- read.csv(paste(dirname(LoaderPATH), "names.csv", sep=.Platform$file.sep))
  for (factor in 2:(dim(keeptable)[2])) { #edited code
    keeprowbool <- sapply(keeptable[,factor], function(x) any(x=="Y"))
    keepnames <- keeptable[keeprowbool,1]
    keeprows <- subset(keeptable, keeprowbool)
    data <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]

    data <- cbind(data, FFdata[,which(colnames(FFdata) %in% colnames(keeptable[-1]))])

    data1 <- apply(data, 2, lgcl) #logicle transform
    message("Applying logicle transform")
    #FFdata1 <- apply(FFdata, 2, lgcl)
    OneDtSNEname <- colnames(keeptable)[factor]
    dataX <- data1[,which (colnames(data1) %in% keeprows[,1])]
    #tSNEdata3 <- Rtsne(dataX, dims=1) #1d rtsne on function markers
    tSNEmat1 <- as.matrix(data[,OneDtSNEname])
    colnames(tSNEmat1) <- OneDtSNEname
    hist(tSNEmat1, 100)
    Xx1DtSNEmat <- cbind(Xx1DtSNEmat,tSNEmat1)

    #Make heatplot annotation


    tSNEbins <- cut(tSNEmat1, breaks = Bins, labels= 1:Bins)  #cut divides the range of x into the number of break intervals specified.
    OneSENSEmed <- apply(dataX, 2, function(x) (tapply(x, tSNEbins, FUN = median))) #tapply(x, index, FUN)
    OneSENSEmed[is.na(OneSENSEmed)] <- 0
    #data2 <- as.data.frame(cbind(data1,tSNEbins))
    #OneSENSEmed <- aggregate(data1 ~ tSNEbins,data =data1, FUN = median, na.rm = F)
    #
    # OneSENSEmed <- data2 %>% group_by(tSNEbins) %>% summarise(median)
    # ddply(data2, "tSNEbins", summarise,
    #       Med = medianmax(year) - min(year),
    #       nteams = length(unique(team)))
    #
    # test <- with(dataX, tapply())


    pdfFN <- paste(dirname(LoaderPATH),paste0(OneDtSNEname,"_Freq.pdf"), sep=.Platform$file.sep)
    pdf(file=pdfFN, width=14, height =6)
    #  my_palette <- colorRampPalette(c("blue","green","yellow"))(n=430)
    breaks = seq(0,3,by=0.005)

    # you can put whatever colors you like here:
    #my_palette <- colorRampPalette(c("black","orange","yellow","white"))(n=length(breaks)-1)
    #my_palette <- colorRampPalette(c("black","blue","green","yellow"))(n=length(breaks)-1)
    #my_palette <- colorRampPalette(c("blue","green","yellow"))(n=length(breaks)-1)
    #my_palette <- colorRampPalette(c("black","blue","green","yellow","darkred"))(n=length(breaks)-1)
    my_palette <- colorRampPalette(c("blue","white","red"))(n=length(breaks)-1)
    #my_palette <- colorRampPalette(c("blue","lightblue","green","yellow","darkorange","darkred"))(n=length(breaks)-1)



    hmap <- heatmap.2(t(OneSENSEmed), col=my_palette,
                      breaks = breaks,
                      margins = c(10,20),
                      Colv = FALSE,
                      dendrogram = "row",
                      cexCol = 1., cexRow =1., scale="none", key=TRUE,trace="none",
                      density.info=c("none"),
                      keysize=1)
    hmapclu= t(OneSENSEmed)[hmap$rowInd, hmap$colInd]
    hmapclut <- t(hmapclu)
    testcol <- seq(from = min(Xx1DtSNEmat), to = (max(Xx1DtSNEmat)-(max(Xx1DtSNEmat)/Bins)), by = (max(Xx1DtSNEmat)/Bins))
    if (OneDtSNEname == "input") {
      p1 <- plot_ly(z= hmapclu,y = rownames(hmapclu), x = testcol,colors = my_palette, type = "heatmap") %>%
        layout(title = paste(OneDtSNEname, "Median heatplot", sep = " "))
      p1
    }
    else {
      p2 <- plot_ly(z= hmapclut,x = colnames(hmapclut), y=testcol,colors = my_palette, type = "heatmap") %>%
        layout(title = paste(OneDtSNEname, "Median heatplot", sep = " "), showlegend = FALSE)
      p2
    }


    dev.off()

  }
  suppressWarnings(OneSplot <<- plot_ly(data.frame(Xx1DtSNEmat), x= Xx1DtSNEmat[, 1], y =Xx1DtSNEmat[, 2], type = "scatter", symbol = "circle-dot"))
  OneSplot
  suppressWarnings(combined <- subplot(OneSplot, p2, p1,nrows = 2, shareY=TRUE, shareX = TRUE))
  combined

  export(combined, file=paste(dirname(LoaderPATH), "groupone.png", sep=.Platform$file.sep))
  browseURL(paste(dirname(LoaderPATH), "groupone.png", sep=.Platform$file.sep))
}


#****************************************************************************************************#
#                                                                                                    #
#****************************************************************************************************#

#' Title
#' Get parameters
#' @param rawFCSdir FCS path
#'
#' @return mid
#' @export
#'
#' @examples
#' #remove hash to run
#' dir3 <- system.file('extdata', package = 'oneSENSE')
#' getParameters(dir3)
getParameters <- function(rawFCSdir) {


  fcsFile <- list.files(path = rawFCSdir, pattern = ".fcs$", full.names = TRUE)
  fcs <- suppressWarnings(read.FCS(fcsFile[1]))
  pd <- fcs@parameters@data
  markers <<- paste("<", pd$name, ">:", pd$desc, sep = "")
  markers_desc <<- pd$desc
  #channels <- pd$name

  if (length(markers) == 0) {
    stop("No markers found in the FCS file!")
  }
}
#****************************************************************************************************#
#                                                                                                    #
#****************************************************************************************************#
#' Title
#'
#' @param LoaderPATH Path of fcs_Out files
#'
#' @return Flow Frame data for coordinate selection
#' @export
#'
#'
#' @examples
#' dir1 <- system.file('extdata',package='oneSENSE')
#' OneSmapperFreq1(dir1) #remove hash symbol to run
OneSmapperFreq1 <- function(LoaderPATH = "fcs") {
  fs <- read.flowSet(path = LoaderPATH, pattern = '.fcs$') #3 sample Fcs files read
  FcsFileNames <- rownames(keyword(fs, "FILENAME")) #store fcs file names to variable fcsfilenames
  FNumBC <- length(fs) #length of fs = 3
  FFdata <- NULL #Flowframedata
  for (FFs in 1:FNumBC) {
    FFa <- exprs(fs[[FFs]]) #Store data from each fcs file to FFa

    #Fixup column names
    colnames(FFa) <- fs[[FFs]]@parameters$desc #assign column names to each column
    empties <- which(is.na(colnames(FFa)) | colnames(FFa)== " ") #find which columns are NA or space
    colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties] #add the corresponding column names to FFa which are initially empty
    fs[[FFs]]@parameters$desc <- colnames(FFa) #add colnames of FFa to the fcs file
    #fs[[FFs]]@parameters$name <- colnames(FFa)
    #Add file label
    FFa <- cbind(FFa,rep(FFs,dim(FFa)[1])) #rep is replicate, it binds the number associated with the file to the last column
    colnames(FFa)[dim(FFa)[2]] <- "InFile" #give the last column the location of the file
    #Concatenate
    FFdata <<- rbind(FFdata,FFa) #adds row data from FFa to FFdata
  }


  lgcl <<- logicleTransform(w=0.05, t=16409, m=4.5, a=0)
  ilgcl <<- inverseLogicleTransform(trans = lgcl)

}

#****************************************************************************************************#
#                                                                                                    #
#****************************************************************************************************#

#' Title
#'
#' @param LoaderPATH Path of fcs_out files
#' @param FFdata Flow frame data
#'
#' @return Assign global variables for access
#' @export
#'
#' @examples
#' #remove hash symbol to run
#' dir4 <- system.file('extdata/extra', package = 'oneSENSE')
#' file5 <- system.file('extdata/FFdata.Rdata', package = 'oneSENSE')
#' load(file5)
#' getCoords(dir4, FFdata=FFdata)
getCoords <- function(LoaderPATH = LoaderPATH, FFdata = FFdata) {
  lgcl <- logicleTransform(w=0.25, t=16409, m=4.5, a=0)
  gckeeptable <- read.csv(paste(dirname(LoaderPATH), "names.csv", sep=.Platform$file.sep))
  gckeeprowbool <- apply(gckeeptable[,c(2,3)],1, function(x) any(x=="Y"))
  gckeepnames <<- gckeeptable[gckeeprowbool,1]

  gckeeprows <- subset(gckeeptable, gckeeprowbool) #return subset of (Function/Diff/Trafficking)
  gcdata <- FFdata[,which (colnames(FFdata) %in% gckeeprows[,1])] # take data that are (Function/Diff/Trafficking)
  gcdata <- cbind(gcdata, FFdata[,which(colnames(FFdata) %in% colnames(gckeeptable[-1]))]) #cbind Function/Diff/Trafficking column from FFdata to data


  gcdata1 <<- apply(gcdata, 2, lgcl) #apply logicle transformation to data


  #Min <- numeric(length = length(keepnames))
  Mid <<- numeric(length = length(gckeepnames))
  #Max <- numeric(length = length(keepnames))
}

#****************************************************************************************************#
#                                                                                                    #
#****************************************************************************************************#

#' Title
#' Frequency Heatplot Generation Part 2
#'
#' @param LoaderPATH path of fcs_Out files
#' @param Bins Number of bins
#'
#' @return returns frequency heatplot
#' @export
#'
#' @examples
#' dir2 <- system.file('extdata',package='oneSENSE')
#'
#' #OneSmapperFreq2(dir2, 250) #remove hash symbol to run
OneSmapperFreq2 <- function(LoaderPATH = "fcs", Bins = 250) {

  Xx1DtSNEmat <- NULL
  keeptable <- read.csv(paste(dirname(LoaderPATH), "names.csv", sep=.Platform$file.sep))
  for (factor in 2:(dim(keeptable)[2])){
    keeprowbool <- sapply(keeptable[,factor], function(x) any(x=="Y"))
    keepnames <- keeptable[keeprowbool,1]
    keeprows <- subset(keeptable, keeprowbool)
    data <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]
    data <- cbind(data, FFdata[,which(colnames(FFdata) %in% colnames(keeptable[-1]))])
    data1 <- apply(data, 2, lgcl) #logicle transform

    OneDtSNEname <- colnames(keeptable)[factor]
    keeprowbool <- sapply(keeptable[,factor], function(x) any(x=="Y"))
    keepnames <- keeptable[keeprowbool,1]
    keeprows <- subset(keeptable, keeprowbool)
    dataX <- data1[,which (colnames(data1) %in% keeprows[,1])]
    tSNEmat1 <- as.matrix(data[,OneDtSNEname])
    colnames(tSNEmat1) <- OneDtSNEname
    hist(tSNEmat1, 100)
    Xx1DtSNEmat <- cbind(Xx1DtSNEmat,tSNEmat1)

    #Make heatplot annotation



    dataGNorm <- dataX
    tSNEbins <- cut(tSNEmat1, breaks = Bins, labels= 1:Bins)   #assigns bin value to each value in tSNEmat1 (1 - 250)

    OneSENSEpp <- matrix(dataX, nrow = Bins, ncol = dim(dataX)[2]) #####@@@@@@ What is the missing data?
    colnames(OneSENSEpp) <- colnames(dataX)
    Coords <- read.csv(paste(dirname(LoaderPATH), "Coords.csv",sep=.Platform$file.sep))
    for(pname in colnames(dataX)) {

      overthresh <- function(group) { #what is this group

        percpos <- (sum(group > Coords[which(Coords$X == pname),2]) / length(group))*100 # calculates percentage of cells that are positive per bin
        #print(percpos)
        return(percpos)
      }

      OneSENSEpp[,pname] <- tapply(dataX[,pname], tSNEbins, FUN = overthresh) #what is tapply
    }

    OneSENSEpp[is.na(OneSENSEpp)] <- 0


    pdfFN <- paste(dirname(LoaderPATH),paste0(OneDtSNEname,"_Freq.pdf"), sep=.Platform$file.sep)
    pdf(file=pdfFN, width=14, height =6)
    #  my_palette <- colorRampPalette(c("blue","green","yellow"))(n=430)
    #breaks = c(seq(0,9,by=0.05),seq(10,100,by=10)) #0.05
    breaks = seq(0,100,by=0.05)

    # you can put whatever colors you like here:
    #my_palette <- colorRampPalette(c("black","orange","yellow","white"))(n=length(breaks)-1)
    #my_palette <- colorRampPalette(c("black","blue","green","yellow"))(n=length(breaks)-1)
    #my_palette <- colorRampPalette(c("blue","green","yellow"))(n=length(breaks)-1)
    #my_palette <- colorRampPalette(c("black","blue","green","yellow","darkred"))(n=length(breaks)-1)
    my_palette <- colorRampPalette(c("blue","white","red"))(n=length(breaks)-1)
    #my_palette <- colorRampPalette(c("blue","lightblue","green","yellow","darkorange","darkred"))(n=length(breaks)-1)



    fhmap <- heatmap.2(t(OneSENSEpp), col=my_palette,
                       breaks = breaks,
                       margins = c(10,20),
                       Colv = FALSE,
                       dendrogram = "row",
                       cexCol = 1., cexRow =1., scale="none", key=TRUE,trace="none",
                       density.info=c("none"),
                       keysize=1)
    fhmapclu= t(OneSENSEpp)[fhmap$rowInd, fhmap$colInd]
    fhmapclut <- t(fhmapclu)
    ftestcol <- seq(from = min(Xx1DtSNEmat), to = (max(Xx1DtSNEmat)-(max(Xx1DtSNEmat)/Bins)), by = (max(Xx1DtSNEmat)/Bins))
    if (OneDtSNEname == "input") {
      suppressWarnings(d1 <- plot_ly(z= fhmapclu,y = rownames(fhmapclu), x = ftestcol,colors = my_palette, type = "heatmap") %>%
                         layout(title = paste(OneDtSNEname, "Frequency heatplot", sep = " ")))
      d1
    }
    else {
      suppressWarnings(d2 <- plot_ly(z= fhmapclut,x = colnames(fhmapclut), y=ftestcol,colors = my_palette, type = "heatmap") %>%
                         layout(title = paste(OneDtSNEname, "Frequency heatplot", sep = " ")))
      d2
    }
    dev.off()

  }
  suppressWarnings(fcombined <- subplot(OneSplot, d2, d1,nrows = 2, shareY=TRUE, shareX = TRUE))
  fcombined
  #tmpFile2 <- tempfile(fileext = ".png")
  export(fcombined, file=paste(dirname(LoaderPATH), "grouptwo.png", sep=.Platform$file.sep))

  browseURL(paste(dirname(LoaderPATH), "grouptwo.png", sep=.Platform$file.sep))

}
