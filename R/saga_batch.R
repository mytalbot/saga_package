#' Batch correction of joint SAGA data sets.
#'
#' \code{saga_batch} removes batch effects from the joint data set. For each sample, the user must specify a batch number in the
#' "SAGA_USER_Samples.txt" file. If no batch numbers are available - you'll still have to provide a number for the samples. You can use
#' any integer number beginning with 1. The batch correction of the data set relies on the ComBat algorithm/package.
#'
#' @param matrix.SAGA.qn quantile normalized saga matrix
#' @param matrix.test.qn quantile normalized test matrix
#' @param rawdata SAGA raw data
#' @param pData.joint matrix with combined SAGA data and user samples.
#' @param plotnumber default=2 for PCA plot with samples in black; alternative: 1 for t-SNE plot
#'
#' @return \code{eset.batch} batch corrected data set.
#'
#' @import sva
#' @import bapred
#' @import Rtsne
#' @importFrom graphics plot abline legend text
#' @importFrom stats prcomp
#'
#' @export
#'


saga_batch     <- function(matrix.SAGA.qn, matrix.test.qn, rawdata, pData.joint, plotnumber=2){

  pData        <- saga::pData
  SIF          <- rawdata$SIF
  pData.Test   <- rawdata$pData.Test
  Top          <- saga::Top11

  ################################################################################################
  #### 3. Addon COMBAT batch correction ##########################################################
  ################################################################################################
  batch.SAGA   <- as.factor(pData$Batch)
  combat.SAGA  <- combatba(t(matrix.SAGA.qn), batch = batch.SAGA)
  matrix.SAGA  <- t(combat.SAGA$xadj)
  colnames(matrix.SAGA) <- row.names(pData)
  matrix.test  <- t(bapred::combatbaaddon(combat.SAGA, t(matrix.test.qn), batch = as.factor(SIF$Batch)))

  #### 3.2 t-SNE of batch corrected dataset ######################################################
  ###################################################################################ä############
  set.seed(12)
  tsne_out     <- Rtsne(t(cbind(matrix.SAGA, matrix.test)), dims = 2, perplexity = 16,
                        theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,
                        verbose = FALSE, is_distance = FALSE)
  #set.seed(59)
  #matrix.joint <- cbind(matrix.SAGA,matrix.test)

  #tsne_out <- Rtsne(t(matrix.joint),dims = 2, perplexity = 16,
  #                  theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,
  #                  verbose = FALSE, is_distance = FALSE)

  if(plotnumber == 2){
    plot(tsne_out$Y,
         col  = pData.joint$IVIM_Farbe,
         pch  = 16,
         cex  = 1.3,
         main = "t-SNE Plot",
         xlab = "Dimension 1",
         ylab = "Dimension 2")
  }else{}


  #### 3.3  clean up data set before assessing SAGA classifier  #################################
  ###############################################################################################
  matrix.train   <- t(matrix.SAGA[row.names(Top),])

  pca.train      <- prcomp(matrix.train, center = T, scale. = T)
  coord.pca.test <- predict(pca.train, newdata = t(matrix.test[row.names(Top),]))

  if(plotnumber == 1){
    plot(rbind(pca.train$x,coord.pca.test),
         pch  = 16,
         col  = pData.joint$Design_Color,
         cex  = 1,
         asp  = 1)
    legend(1,-2, legend = c("transforming","mock","neutral","new samples"), col = unique(pData.joint$Design_Color), pch=16, bty="n", cex=0.8)
    text(coord.pca.test[,c(1:2)], labels=pData.Test$Filename, cex= 0.3, pos=3, offset = 0.3)
  }else{}

  index        <- nrow(pData)+nrow(pData.Test)


  ##matrix.Top   <- cbind(matrix.SAGA, matrix.test)[row.names(Top),]
  ##index        <- nrow(pData)+nrow(pData.Test)


  ##if(plotnumber == 1){
  ##  pca     <- prcomp(t(matrix.Top))
  ##  plot(pca$x,
  ##       pch  = 16,
  ##       col  = pData.joint$Design_Color,
  ##       cex  = 1,
  ##       asp  = 1,
  ##       xlab = "Dimension 1",
  ##       ylab = "Dimension 2")
  ##  legend("topleft", legend = c("transforming","mock","neutral","new samples"), col = unique(pData.joint$Design_Color), pch=16, bty="n", cex=0.8)
  ##  text(pca$x[c((nrow(pData)+1):index),c(1:2)] , labels=pData.Test$Filename, cex= 0.6, pos=3, offset = 0.3) # new in V6: (nrow(pData)+1)
  ## }else{}




  return(list(matrix.SAGA = matrix.SAGA, matrix.test = matrix.test, index = index))

}


