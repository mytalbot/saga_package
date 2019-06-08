#' Prediction of sample arrays with the SAGA model.
#'
#' \code{saga_predict} uses a Support Vector Machine (SVM) with a radial kernel to classify user samples either as transforming or nontransforming.
#' The SVM model is built with the integrated SAGA training data set and a toplist of relevant probes which are able to differentiate assay data
#' into the mentioned risk states. The SVM is further optimized on the SAGA training data. Use at your own risk!
#'
#' @param samplepath path to sample data
#' @param matrix.train normalized, probe-averaged and batch-corrected SAGA training data.
#' @param labels.train class labels (factors) for SAGA training data. Can either be "transforming" or "nontransforming".
#' @param matrix.unknown matrix of sample data with array names as row names and probes as column names.
#' @param pData.Test SAGA sample or test data matrix
#' @param writeFile default=1 for writing results to the sample folder
#' @param showRoc default=1 for showing the ROC curve of the model performance; alternative: 0 for showing naught
#'
#' @return \code{predictions} Data frame with three columns. Column one shows the sample names, the second column shows the decision values of
#' the svm function and column thress shows the predictions for the query assays either as transforming or nontransforming.
#'
#' @import caret
#' @import pROC
#' @import lattice
#' @import ggplot2
#' @import utils
#' @importFrom stats predict
#'
#' @export
#'


saga_predict <- function(samplepath, matrix.train, labels.train, matrix.unknown, pData.Test, writeFile=0, showRoc=1) {

  ################################################################################################
  #### 6. Caret SVM ##############################################################################
  ################################################################################################
  oldw <- getOption("warn")
  options(warn = -1)

  #### 6.1 Caret SVM on ClassProbs ###############################################################
  ################################################################################################
  set.seed(45)
  svm_fit.ROC <- train(matrix.train,
                       labels.train,
                       method     = "svmRadial",
                       tuneLength = 20,
                       trControl  = trainControl(method     = "repeatedcv",number  = 10, repeats = 5,
                                                 classProbs = TRUE,summaryFunction = twoClassSummary))
  print(svm_fit.ROC)

  Prediction_SVM.Caret <- predict(svm_fit.ROC, matrix.unknown, type = "prob")
  Prediction_SVM.Caret$Prediction.SVM.Caret <- ifelse(Prediction_SVM.Caret$transforming>0.50,"transforming","untransforming")
  Prediction_SVM.Caret <- cbind(pData.Test[,c(1:5)],Prediction_SVM.Caret)


  if(writeFile==TRUE){
    write.table(Prediction_SVM.Caret, file = paste(samplepath,"/result_SAGA_Predictions.txt",sep = ""),
                sep="\t", row.names = TRUE, col.names=NA, quote=FALSE)

  #### 6.2 Performance of Classifier for n samples ##############################################
  ################################################################################################
  sink( paste(samplepath,"/result_ConfusionMatrix_samples.txt", sep = ""), append = TRUE)
    print( confusionMatrix(as.factor(Prediction_SVM.Caret$Prediction.SVM.Caret), as.factor(Prediction_SVM.Caret$TrueLabel)) )
  sink()

  }else{}

  if(showRoc==1){
    Prediction_SVM.Caret$Class <- as.factor(ifelse(Prediction_SVM.Caret$TrueLabel == "transforming","transforming","nontransforming"))

    roc1 <- roc(Prediction_SVM.Caret$Class,                    # response vector (factor or character)
                Prediction_SVM.Caret$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T,print.thres=F, main="ROC SAGA Samples")
  }else{}


  options(warn = oldw)

  return(list(svm_fit=svm_fit.ROC, predictions=Prediction_SVM.Caret))

}








