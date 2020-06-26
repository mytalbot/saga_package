#' Geneset Enrichment Analysis (GSEA) of SAGA data.
#'
#' \code{saga_gesea} performs batch-wise GSEA for SAGA samples.
#'
#' @param smplpath path to the saga data folder with the user samples.
#' @param saveResults can be 1 for yes and 0 for no.
#'
#' @return \code{result} GSEA result - will also be saved into the sample folder.
#'
#' @import limma
#' @importFrom gridExtra combine grid.table
#' @importFrom stats setNames
#' @importFrom methods new
#' @importFrom graphics plot abline legend text
#' @import utils
#' @importFrom grDevices pdf dev.off
#' @import phenoTest
#' @importFrom Biobase ExpressionSet
#' @import GSEABase
#'
#' @export
#'

saga_gesea    <- function(smplpath, saveResults=0){

  sets        <- saga::sets
  SAGA.CORE   <- setNames(split(sets, seq(nrow(sets))), rownames(sets))   # GeneSets have to be stored in a list object
  SIF         <- read.delim(paste(smplpath,"/SampleInformation.txt",sep=""),
                            row.names=1,header=TRUE,sep="\t", stringsAsFactors = F) # Sample Information File for phenoTest


  ### 1.1. Read in files from test and loop over all batches separately from here on ##############
  #################################################################################################
  maxBatch   <- max(as.integer(SIF$Batch))   # how many assays / batches


  for(i in 1:maxBatch) {
    SIF.i <- SIF[SIF$Batch==i,]
    RAW.i <- read.maimages(files=SIF.i$Filename, path=smplpath, source="agilent.median", green.only=T,
                           columns=list(G="gMedianSignal"), annotation=c("ProbeName", "GeneName"))
    colnames(RAW.i) <- row.names(SIF.i)
    #### 2.1. Normalize, average ###################################################################
    RMA.i <- normalizeBetweenArrays(RAW.i, method="quantile")            # quantil normalize
    RMA.i <- avereps(RMA.i,ID= RMA.i$genes$ProbeName)                    # average replicates to one value for each probe
    matrix.gsea <- RMA.i$E                                               # extract log2 expression values

    #### 2.3. make ExpressionSet (Biobase) object ##################################################
    metadata  <- data.frame(labelDescription= rep(NA,dim(SIF.i)[2]),row.names=colnames(SIF.i))   # varMetadata: empty, but required
    phenoData <- new("AnnotatedDataFrame",data=SIF.i, varMetadata = metadata)     # annotatedDataFrame for the annotation of the samples
    eset.gsea <- ExpressionSet(assayData = matrix.gsea, phenoData = phenoData)  # this is the ExpressionSet required for phenoTest

    #### 2.4. make ePheno object: contains the FCs associated with Group variable ##################
    vars2test   <- list(ordinal="Group")    # Variables (here: Groups) to test against MOCK, which are always Group = 1 in the SIF
    epheno.gsea <- ExpressionPhenoTest(eset.gsea,vars2test,p.adjust.method='BH')

    #### 2.5 GSEA #################################################################################
    SAGA.GSEA <- gsea(x=epheno.gsea, gsets=SAGA.CORE ,B=2000,                  # calculate GSEA-scores based on the FC in the epheno object
                      center = TRUE, test = "perm", p.adjust.method='BH')

    result            <- summary(SAGA.GSEA)[,c(1,2,3,5,8)]                     # extract results (only NES- normalized enrichment scores)
    #NEW
    result$pred.class <- ifelse(result$nes>0.9,"transforming","untransforming")  # prediction based on NES

    #### 2.6 output ###############################################################################
    Group <- NULL    ### pull out the Group index number from the result table
    for (a in 1:nrow(result)) {Group[a] <- unlist(strsplit(as.character(result$variable[a]), ".", fixed = TRUE))[2] }
    result$Group      <- Group

    SIF.sub           <- SIF.i[SIF.i$Group != 1, c(3,4,1) ]                  # pull out info of tested Groups
    SIF.sub$SampleID  <- row.names(SIF.sub)
    result.m          <- merge(SIF.sub,result, by.x="Group", by.y = "Group") # merge result with SIF for SampleIDs and FileNames
    write.table(result.m, file = paste(smplpath,"/Results_SAGA.GSEA_Batch_",i,".txt",sep = ""), sep="\t",row.names = FALSE)

    # make pdf report
    pdf(file=paste(smplpath,"/","SAGA.GSEA_Batch_",i,".pdf",sep = ""),useDingbats = F,width = 10, height = 10)
    gridExtra::grid.table(result.m,rows = NULL)
    plot(SAGA.GSEA,es.nes='nes',selGsets='SAGA.CORE')
    dev.off()
  }
}




