# Evan Newell 2015
# Edited by Tan Yong Kee

#' tSNE and OneSENSE algorithm for FCS data
#'
#' @param LoaderPATH Path where FCS file is located
#' @param ceil Maximum number of cells to sample from each fcs sample/file
#' @param FNnames .csv file generated when markers from each category are selected
#' @param OutputSuffix suffix to name output folder
#' @param DotSNE boolean, if TRUE do tSNE, if FALSE skip tSNE
#' @param DoOneSENSE boolean, if TRUE do OneSENSE, if FALSE skip OneSENSE
#' @param Bins number of bins to put the cell data into, DEFAULT = 250
#'
#' @return FCS files, tSNE histograms, OneSENSE plot
#'
#' @importFrom Rtsne Rtsne
#' @importFrom graphics hist
#' @importFrom flowCore read.flowSet exprs keyword write.FCS logicleTransform inverseLogicleTransform identifier exprs<- identifier<-
#' @importFrom methods cbind2
#'
#' @examples
#' #dir <- system.file('extdata',package='oneSENSE')
#' #fnnames <- system.file('extdata', "names.csv", package='oneSENSE')
#'
#' #FCStSNE(LoaderPATH = dir, FNnames = fnnames) #remove hash symbol to run
FCStSNE <- function (LoaderPATH ="fcs",
                     ceil = 5000,
                     FNnames="names.csv",
                     OutputSuffix = "Out",
                     DotSNE = TRUE,
                     DoOneSENSE = TRUE,
                     Bins = 250) #it doesnt use the variable bins in the code below
{

  fs <-read.flowSet(path = LoaderPATH, pattern = ".fcs$") #Read one or several FCS files, path is the directory
  FcsFileNames <- rownames(keyword(fs, "FILENAME")) #keyword: Retreive keywords of a flowFrame #rownames require
  NumBC <- length(fs) #3
  FFdata <- NULL #FlowFrame data

  for (FFs in 1:NumBC){ #iterate through each FCS file
    FFt <- exprs(fs[[FFs]]) # assign matrix (row:events, columns (55 or them) are parameters (time, event length, metal))
    ## Downsample ##
    if (nrow(FFt)<=ceil) #number of events less than assigned ceil, assign to FFa
      FFa <- FFt
    else
      FFa <- FFt[sample(nrow(FFt),ceil,replace=FALSE),] #sample takes a sample of the specified size form the elements of x using either with or without replacement
    #in this case, a specified number of rows are taken randomly from the matrix without replacement
    #Fixup column names
    colnames(FFa) <- fs[[FFs]]@parameters$desc # set the column names of the matrix by retrieving the parameters (protein markers)
    empties <- which(is.na(colnames(FFa)) | colnames(FFa)== " ") # assign empties to the empty columns (type: integer)
    colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties] # write in col names that are NA with appropriate column names (Time, Event_length, Length)
    fs[[FFs]]@parameters$desc <- colnames(FFa) #assign column names @parameters$desc from FFa to fs[[1,2,3]]
    #fs[[FFs]]@parameters$name <- colnames(FFa)
    #Add file label
    FFa <- cbind(FFa,rep(FFs,dim(FFa)[1])) #add a new column with each row having an entry indicating the sample (1,2,3)
    colnames(FFa)[dim(FFa)[2]] <- "InFile" #name the 56th column "InFile"
    #Concatenate
    FFdata <- rbind(FFdata,FFa) #what's the difference between FFdata and FFa, seems the same to me, is it a clone? Has all the markers
    #rbind concatenates the rows to its respective rows. All 3 samples of 20000 rows together to make 60000 rows
  }
  message("FCS Files Read")
  keeptable <- read.csv(paste(dirname(LoaderPATH), "names.csv", sep=.Platform$file.sep)) #assign read csv file to keeptable
  keeprowbool <- sapply(keeptable[,2], function(x) any(x=="Y")) #apply a function x to each element of a vector, return a row of booleans(if y is present its true, else false)
  keeprows <- subset(keeptable, keeprowbool) #Return subsets of vectors, matrices or data frames which meet conditions. IN this case only function markers
  data <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])] #new data file with only function markers

  #lgcl <- logicleTransform(w=0.1, t=500000, m=4.5, a=0) #hyperbolic sine transformation functions for display of flow cytometry data.
  lgcl <- logicleTransform(w=0.25, t=16409, m=4.5, a=0)
  ilgcl <- inverseLogicleTransform(trans = lgcl)

  data1 <- apply(data, 2, lgcl) #apply logicle transform to data1 (function markers), 2 means apply over columns
  FFdata1 <- apply(FFdata, 2, lgcl) #apply to all markers


  score <- NULL
  tSNEmat <- NULL
  if(DotSNE) # does tSNE for the 15 functional markers only.
  {
    message("Doing tSNE")
    tSNEdata3 <- Rtsne(data1, dims=2) #returns a list of dimension reduced data
    tSNEmat <- tSNEdata3$Y #Matrix containing the new representations for the objects (all)
    colnames(tSNEmat) <- c("tSNE1","tSNE2") #assignment col names to tSNEmat
    plot(tSNEmat[, 1], tSNEmat[, 2], pch=".", xlab="tSNE1", ylab="tSNE2", cex=0.1)
  } else tSNEmat <- NULL


  if(DoOneSENSE)
  {
    message("Doing oneSENSE")
    Xx1DtSNEmat <- NULL
    if (dim(keeptable)[2] == 4) {
      for (factor in 2:(dim(keeptable)[2])){ #for loop from 2 to 4
        OneDtSNEname <- colnames(keeptable)[factor]
        keeprowbool <- sapply(keeptable[,factor], function(x) any(x=="Y"))
        keeprows <- subset(keeptable, keeprowbool)
        dataX <- FFdata1[,which (colnames(FFdata1) %in% keeprows[,1])] #apply one dimension tsne to each catergory of cells per loop
        tSNEdata3 <- Rtsne(dataX, dims=1, check_duplicates=FALSE)
        tSNEmat1 <- tSNEdata3$Y
        colnames(tSNEmat1) <- OneDtSNEname
        hist(tSNEmat1, 100)
        Xx1DtSNEmat <- cbind(Xx1DtSNEmat,tSNEmat1)
      }
      #plot(Xx1DtSNEmat[, 1], Xx1DtSNEmat[, 2], pch=".", xlab=paste("tSNE1 ",colnames(Xx1DtSNEmat)[1],sep=""), ylab=paste("tSNE2 ",colnames(Xx1DtSNEmat)[2],sep=""), cex=1)
      scatterplot3d(x = Xx1DtSNEmat[, 1], y = Xx1DtSNEmat[, 2], z= Xx1DtSNEmat[,3], pch=".",
                    xlab=paste("tSNE1 ",colnames(Xx1DtSNEmat)[1],sep=""),
                    ylab=paste("tSNE2 ",colnames(Xx1DtSNEmat)[2],sep=""),
                    zlab=paste("tSNE3 ",colnames(Xx1DtSNEmat)[3],sep="")
      )
    }
    else { #dim(keeptable)[2] == 3
      for (factor in 2:(dim(keeptable)[2])){ #for loop from 2 to 4
        OneDtSNEname <- colnames(keeptable)[factor]
        keeprowbool <- sapply(keeptable[,factor], function(x) any(x=="Y"))
        keeprows <- subset(keeptable, keeprowbool)
        dataX <- FFdata1[,which (colnames(FFdata1) %in% keeprows[,1])] #apply one dimension tsne to each catergory of cells per loop
        tSNEdata3 <- Rtsne(dataX, dims=1, check_duplicates=FALSE)
        tSNEmat1 <- tSNEdata3$Y
        colnames(tSNEmat1) <- OneDtSNEname
        hist(tSNEmat1, 100)
        Xx1DtSNEmat <- cbind(Xx1DtSNEmat,tSNEmat1)
      }
      plot(Xx1DtSNEmat[, 1], Xx1DtSNEmat[, 2], pch=".",
           xlab=paste("tSNE1 ",colnames(Xx1DtSNEmat)[1],sep=""),
           ylab=paste("tSNE2 ",colnames(Xx1DtSNEmat)[2],sep=""),
           cex=1)
    }
  }
  else Xx1DtSNEmat <- NULL


  NXx1 <- apply(Xx1DtSNEmat,2,function(x) ((x-min(x))/(max(x)-min(x)))*10000 ) #consider each column, normalise each value in the column.

  score2 <- cbind(score,tSNEmat)
  Nscore <- apply(score2,2,function(x) ((x-min(x))/(max(x)-min(x)))*3.7 ) #normalise tSNEmat scores
  ilgcl <- inverseLogicleTransform(trans = lgcl)
  NIscore <- apply(Nscore, 2, ilgcl)
  colnames(NIscore) <- colnames(score2)

  NIscore <- cbind(NIscore, NXx1)

  message("Writing FCS Files")
  #output new FCS files
  for (FFs in 1:NumBC){

    newFF <- fs[[1]]
    newBaseData <- FFdata[FFdata[,dim(FFdata)[2]]==FFs,-dim(FFdata)[2]]
    colnames(newBaseData)<- colnames(exprs(newFF))
    exprs(newFF) <- newBaseData
    subsetNIscore <- NIscore[FFdata[,dim(FFdata)[2]]==FFs,]
    newFF <- cbind2(newFF, subsetNIscore)
    newFF@parameters$desc <- colnames(cbind(newBaseData, subsetNIscore))
    suppressWarnings(dir.create(paste0(LoaderPATH,"_",OutputSuffix)))
    BaseFN <- sapply(strsplit(FcsFileNames[FFs], split ="\\."), "[", 1)
    FNresult <- paste0(LoaderPATH,"_",OutputSuffix,"/",BaseFN,"_",OutputSuffix,".fcs")
    newFF@description$'$FIL' <- paste0(BaseFN,"_",OutputSuffix,".fcs")
    newFF@description$FILENAME <- paste0(BaseFN,"_",OutputSuffix,".fcs")
    identifier(newFF) <- paste0(BaseFN,"_",OutputSuffix)


    suppressWarnings(write.FCS(newFF, FNresult))
  }
  message("New FCS Files Written")
}
