
#' Sample preparation of the joint SAGA data set.
#'
#' \code{saga_sampling} will prepare training and testing data of the joint SAGA data set. An internal process
#' is used to filter all data for candidate probes suited for classification of array data into transforming/untransforming.
#' The list cannot be changed in this version and the package is optimized for this specific gene set. Also, a PCA plot may be
#' generated to monitor the generalized context of the sample data within the SAGA training set.
#'
#' @param matrix.SAGA Matrix with SAGA model data
#' @param matrix.test Matrix with user sample or test data
#'
#' @return \code{matrix.train} normalized, probe-averaged and batch-corrected SAGA training data.
#' @return \code{labels.train} class labels (factors) for SAGA training data. Can either be "transforming" or "untransforming".
#' @return \code{matrix.unknown} matrix of sample data with array names as row names and probes as column names.
#'
#' @export
#'

saga_sampling    <- function(matrix.SAGA, matrix.test ){


  ################################################################################################
  #### 5. split into Prediction and Known Sets ###################################################
  ################################################################################################
  Top            <- saga::Top12
  pData          <- saga::pData

  matrix.train   <- t(matrix.SAGA[row.names(Top),])
  labels.train   <- as.factor(pData$TrueLabel)
  matrix.unknown <- t(matrix.test[row.names(Top),])


   return(list(matrix.train   = matrix.train,
               labels.train   = labels.train,
               matrix.unknown = matrix.unknown ))

}
