#' Import user data for SAGA model building
#'
#' The \code{saga_import} function loads SAGA data from textfiles in a user specified path and prepares the necessary objects for later
#' classification of unknown arrays.
#'
#' In addition to the sample files, a user-defined targets file must be provided. It contains the necessary sample names, batch
#' information etc. You can use the \code{saga_targets} function to create a targets file for you. However, in any case, the
#' targets file must be named "SAGA_USER_Samples.txt" and must be placed in the same folder together with the sample files.
#'
#' \cr PLEASE NOTE: The saga package works with raw Agilent data. There is no need for prior data preprocessing.
#' The saga function will take care of this.
#'
#' @param smplpath path to the saga data folder with the user samples.
#' @param showjoint can be 1 or 0. In case of 1 it will show a boxplot of the unnormalised data (model + user data). Default is 0.
#'
#' @return \code{matrix.joint} n-by-m matrix of joint data (SAGA data + user samples).
#' @return \code{pData.joint} joint target matrix of both SAGA data and user samples.
#' @return \code{pData.user} target matrix of user samples.
#' @return \code{filelist} list of user array names.
#' @return \code{pData} phenotype target matrix of SAGA data
#'
#' @import limma
#'
#' @export
#'

saga_import   <- function(smplpath, showjoint=0){


  # catch some errors
  #path                <- getwd()

  ################################################################################################
  #### 1. SAGA data: 91 training samples and 9 predictors found by genetic algorithm #############
  ################################################################################################
  # load internal data
  Annotation          <- saga::Annotation
  pData               <- saga::pData
  #targets             <- saga::saga_targets

  SAGA_Data           <- read.delim(system.file("extdata", "SAGA_INBUILD_Data_AVE_91.txt", package = "saga"),header=TRUE,sep="\t", stringsAsFactors =FALSE)
  SAGA_RAW            <- as.matrix(SAGA_Data[,-1] )
  row.names(SAGA_RAW) <- SAGA_Data$PROBE_ID


  ### Read in .txt files from validation set: ####################################################
  ################################################################################################
  SIF                     <- read.delim(paste(smplpath,"/SampleInformation.txt",sep=""),row.names=1,header=TRUE,sep="\t", stringsAsFactors = F) # Sample Information File for phenoTest
  pData.Test              <- read.delim(paste(smplpath,"/SampleInformation.txt",sep=""),row.names=1,header=TRUE,sep="\t", stringsAsFactors = F) # Sample Information File for SAGA-SVM
  pData.Test$Batch        <- pData.Test$Batch + 11              # 11 Batches are already in the SAGA Inbuild DataSet
  pData.Test$IVIM_Color   <- rep("#000000", nrow(pData.Test))
  pData.Test$Design_Color <- rep("#000000", nrow(pData.Test))
  pData.Test$IVIM_ID      <- rep(NA, nrow(pData.Test))
  pData.Test$IVIM_Farbe   <- rep("#000000", nrow(pData.Test))

  TEST_Data               <- read.maimages(files=pData.Test$Filename, path=smplpath, source="agilent.median", green.only=T,
                                           columns=list(G="gMedianSignal"), annotation=c("ProbeName", "GeneName"))

  colnames(TEST_Data)     <- row.names(pData.Test)
  TEST_RAW                <- TEST_Data$E                                       # export Expression matrix
  row.names(TEST_RAW)     <- TEST_Data$genes$ProbeName                         # assign PROBE IDs as row.names
  TEST_RAW                <- avereps(TEST_RAW, ID= row.names(TEST_RAW))        # collapse quadruplicate probes
  TEST_RAW                <- TEST_RAW[row.names(SAGA_RAW),]                    # Expressionmatrix of Test Set with 36226 annotated probes


  ### Make joint sample information file #########################################################
  ################################################################################################
  all( colnames(pData.Test) == colnames(pData) )
  pData.joint <- rbind(pData,pData.Test)

  matrix.joint<- cbind(SAGA_RAW,TEST_RAW)

  ### test samples are in black ##
  if(showjoint ==1)
    boxplot(log2(cbind(SAGA_RAW,TEST_RAW)),
            col       = pData.joint$IVIM_Farbe,
            names     = pData.joint$Filename,
            boxwex    = 0.6,
            cex.axis  = 0.5,
            las       = 2,
            outline   = FALSE,
            main      = "SAGA joint data set (raw)")

  return( list(SAGA_RAW= SAGA_RAW, TEST_RAW=TEST_RAW, pData.joint= pData.joint, pData.Test= pData.Test, pData= pData, SIF=SIF) )

}



