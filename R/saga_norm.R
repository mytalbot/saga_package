#' Normalize joint SAGA data sets
#'
#'
#' \code{saga_norm} uses both, user data and SAGA data to perform quantile normalization. Multiple probes within arrays are averaged.
#'
#' @param SAGA_RAW SAGA raw data matrix
#' @param TEST_RAW Test or sample data in SAGA matrix format
#' @param pData.joint joint data matrix (internal SAGA data + sample data)
#' @param plotnumber default=2 for t-SNE of SAGA model data + sample; alternative = 1 for quantile normalized boxplot with samples
#'
#' @return \code{eset.rma} expression set of normalized array data.
#'
#' @import limma
#' @import UsingR
#' @import bapred
#' @import Rtsne
#'
#' @export
#'

saga_norm <- function(SAGA_RAW, TEST_RAW, pData.joint, plotnumber=2) {

  ################################################################################################
  #### 2.  Addon Quantile Normalization   ########################################################
  ###################################################################################ä############
  qunorm.SAGA     <- qunormtrain(t(SAGA_RAW))
  matrix.SAGA.qn  <- log2(t(qunorm.SAGA$xnorm))
  matrix.test.qn  <- log2(t(qunormaddon(qunorm.SAGA, t(TEST_RAW))))

  ### test samples are in black ##
  if(plotnumber == 1){
   boxplot(cbind(matrix.SAGA.qn, matrix.test.qn),
          boxwex   = 0.6,
          col      = pData.joint$IVIM_Farbe,
          names    = pData.joint$Filename,
          cex.axis = 0.5,
          las      = 2,
          outline  = FALSE,
          main     = "SAGA joint data set (normalized)")
  }else{}


  matrix.joint.qn <- cbind(matrix.SAGA.qn,matrix.test.qn)

  ## for robustness of t-SNE prior dimensionality reduction with PCA is performed (initial dimensions = 50 (default))
  set.seed(12)
  tsne_out.qn <- Rtsne(t(matrix.joint.qn),dims = 2, perplexity = 16,
                       theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,
                       verbose = FALSE, is_distance = FALSE)

  if(plotnumber == 2){  ## test samples are in black
    plot(tsne_out.qn$Y,
         col  = pData.joint$IVIM_Farbe,
         xlim = c(-14,14),
         ylim = c(-14,14),
         pch  = 16,
         cex  = 1.3,
         xlab = "Dimension 1",
         ylab = "Dimension 2")

    text(tsne_out.qn$Y, pData.joint$Filename, cex=0.4, offset = 0.3, pos=3)
  }else{}


return(list(qunorm.SAGA    = qunorm.SAGA,
            matrix.SAGA.qn = matrix.SAGA.qn,
            matrix.test.qn = matrix.test.qn))

}








