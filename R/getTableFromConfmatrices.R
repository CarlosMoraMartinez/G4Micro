#' @title Extract Classification Metrics from Model List
#' @description Takes a list of trained classification models and extracts key evaluation metrics
#' (e.g., Accuracy, Kappa, AUC, Precision, Recall, etc.) from their confusion matrices and ROC AUC results.
#' It summarizes them in a tidy data frame for easy comparison.
#'
#' @param modlist A named list of models as returned by the training and evaluation pipeline. Each model must contain the confusion matrices (`confmat` for cross-validated results and `confmat_no_l1o` for non-cross-validated results) and optionally the ROC AUC (`roc_auc`).
#'
#' @return A data frame containing performance metrics (Accuracy, Kappa, AUC, Sensitivity, Specificity, etc.) for each model, both with and without cross-validation.
#'
#' @details This function extracts evaluation metrics from a list of classification model objects and returns a summary table. Each model in the list must include confusion matrix results (both with cross-validation and without), as well as ROC AUC values, if available. This is typically used after performing K-fold cross-validation (not limited to leave-one-out) on classification models. The function requires that each element of the list is named and that all model results are structured consistently.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}
#' @rdname getTableFromConfmatrices
#' @export
#' @importFrom dplyr mutate select
#' @importFrom caret confusionMatrix
getTableFromConfmatrices <- function(modlist){
  res <- map(modlist, \(mod){
    data.frame(
      Accuracy_l1out=if(is.null(mod$confmat)) NA else mod$confmat$overall["Accuracy"],
      Kappa_l1out=if(is.null(mod$confmat)) NA else mod$confmat$overall["Kappa"],
      Sensitivity_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Sensitivity"],
      Specificity_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Specificity"],
      PPV_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Pos Pred Value"],
      NPV_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Neg Pred Value"],
      Precision_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Precision"],
      Recall_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Recall"],
      BalancedAccuracy_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Balanced Accuracy"],
      AUC_l1out = if(is.null(mod$roc_auc)) NA else mod$roc_auc,
      Accuracy=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Accuracy"],
      Kappa=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Kappa"],
      Sensitivity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Sensitivity"],
      Specificity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Specificity"],
      PPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Pos Pred Value"],
      NPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Neg Pred Value"],
      Precision = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Precision"],
      Recall = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Recall"],
      BalancedAccuracy = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Balanced Accuracy"]
    )
  }) %>% bind_rows() %>%
    dplyr::mutate(model = names(modlist)) %>%
    dplyr::select(model, everything()) %>%
    arrange(desc(Accuracy_l1out))
  rownames(res) <- NULL
  return(res)
}
