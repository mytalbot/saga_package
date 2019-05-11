#' SAGA Wrapper
#'
#' The \code{saga_wrapper} function integrates all SAGA functions into one. With this function you are not as flexible as
#' with the single ones but you can run the analysis/ classification in one step.
#' However, this requires a complete SampleInformatin.txt file. For GESEA analysis the SampleInformation file must also fulfill
#' the neccessary criteria. See saga_vignette for more information.
#'
#' \code{saga_wrapper}
#'
#' @param samplepath sample path
#' @param doGESEA Can be 0 or 1. If set to 1, GESEA analysis will be added to the analysis.
#' @param showModel Can be 0 or 1. If set to 1 (default) PCA plot (+samples) will be shown.
#'
#' @return \code{output} optimized SVM model for classification.
#'
#' @import limma
#' @import bapred
#' @import sva
#' @import UsingR
#' @import caret
#' @import Biobase
#' @import phenoTest
#' @import gridExtra
#'
#' @export
#'

saga_wrapper   <- function(samplepath, showModel=1, doGESEA=0){

  # saga_import
  rawdata        <- saga_import(samplepath, showjoint=0)

  SAGA_RAW       <- rawdata$SAGA_RAW
  TEST_RAW       <- rawdata$TEST_RAW
  pData.joint    <- rawdata$pData.joint
  pData.Test     <- rawdata$pData.Test
  print("+++ sample data read +++")

  # normalize saga data (plotnumber=1 for normalized boxplot, 2 for tSNE plot)
  normalized     <- saga_norm(SAGA_RAW, TEST_RAW, pData.joint, plotnumber=1)

  qunorm.SAGA    <- normalized$qunorm.SAGA
  matrix.SAGA.qn <- normalized$matrix.SAGA.qn
  matrix.test.qn <- normalized$matrix.test.qn
  print("+++ normalization complete +++")

  # remove batch effects
  if(showModel==1){
    batchnorm    <- saga_batch(matrix.SAGA.qn, matrix.test.qn, rawdata, pData.joint, plotnumber=1)
  }else{
    batchnorm    <- saga_batch(matrix.SAGA.qn, matrix.test.qn, rawdata, pData.joint, plotnumber=2)
  }

  index          <- batchnorm$index
  matrix.SAGA    <- batchnorm$matrix.SAGA
  matrix.test    <- batchnorm$matrix.test
  print("+++ batch effects removed +++")

  # Sampling and model data collection
  model          <- saga_sampling(matrix.SAGA, matrix.test )

  matrix.train   <- model$matrix.train
  labels.train   <- model$labels.train
  matrix.unknown <- model$matrix.unknown

  # Array predictions with optimized SVM parameters (default settings)
  output         <- saga_predict(samplepath, matrix.train, labels.train, matrix.unknown,
                                 pData.Test, writeFile=1, showRoc=0)
  print("+++ SVM prediction complete +++")
  print("+++ Prediction files (results) written +++")

  # Predictions with own neg/pos controls and grid optimization; may be used without controls
  classes        <- output$predictions
  print(output)

  print("+++ SVM grid optimization & prediction complete +++")
  print("+++ SAGA ANALYSIS COMPLETE! +++")

  if(doGESEA == 1){
    gesea_results  <- saga_gesea(samplepath, saveResults=1)
    print("+++ GESEA COMPLETE! +++")
    return(list(predictions=output, gesea=gesea_results ) )
  }else{
    return(list(predictions=output ) )
  }

}
