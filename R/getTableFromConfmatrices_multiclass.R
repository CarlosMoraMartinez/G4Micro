#' @title Summarize Performance Metrics from Multiclass Classification Models
#' @description akes a list of trained classification models and extracts key evaluation metrics
#' (e.g., Accuracy, Kappa, AUC, Precision, Recall, etc.) from their confusion matrices and ROC AUC results.
#' It summarizes them in a tidy data frame for easy comparison.
#'
#' @param modlist A named list of trained model objects. Each element should contain confusion matrices
#' (for cross-validated and full model) and optionally a `roc_auc` value. Expected fields are:
#' \itemize{
#'   \item \code{confmat}: Confusion matrix object (e.g., from caret::confusionMatrix) for cross-validation folds
#'   \item \code{confmat_no_l1o}: Confusion matrix for the full training predictions (non-cross-validated)
#'   \item \code{roc_auc}: Optional AUC score for the cross-validation evaluation
#' }
#'
#' @return A data frame summarizing the mean values of multiple classification performance metrics
#' including Accuracy, Kappa, Sensitivity, Specificity, PPV, NPV, Precision, Recall, Balanced Accuracy, and AUC.
#' Metrics are provided separately for cross-validation (\code{*_l1out}) and for full data (\code{*}).
#'
#' @details This function is designed for evaluating multiclass classifiers trained using general cross-validation.
#' It averages the per-class metrics (e.g., Sensitivity, Specificity) across all classes in a consistent way.
#'
#'@examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}
#' @rdname getTableFromConfmatrices_multiclass
#' @export
#' @importFrom dplyr mutate select
getTableFromConfmatrices_multiclass <- function(modlist){
  res <- map(modlist, \(mod){
    data.frame(
      Accuracy_l1out=if(is.null(mod$confmat)) NA else mod$confmat$overall["Accuracy"] %>% mean,
      Kappa_l1out=if(is.null(mod$confmat)) NA else mod$confmat$overall["Kappa"] %>% mean,
      Sensitivity_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Sensitivity"] %>% mean,
      Specificity_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Specificity"] %>% mean,
      PPV_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Pos Pred Value"] %>% mean,
      NPV_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Neg Pred Value"] %>% mean,
      Precision_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Precision"] %>% mean,
      Recall_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Recall"] %>% mean,
      BalancedAccuracy_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Balanced Accuracy"] %>% mean,
      AUC_l1out = if(is.null(mod$roc_auc)) NA else mod$roc_auc,
      Accuracy=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Accuracy"] %>% mean,
      Kappa=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Kappa"] %>% mean,
      Sensitivity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Sensitivity"] %>% mean,
      Specificity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Specificity"] %>% mean,
      PPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Pos Pred Value"] %>% mean,
      NPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Neg Pred Value"] %>% mean,
      Precision = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Precision"] %>% mean,
      Recall = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Recall"] %>% mean,
      BalancedAccuracy = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat$byClass[, "Balanced Accuracy"] %>% mean
    )
  }) %>% bind_rows() %>%
    dplyr::mutate(model = names(modlist)) %>%
    dplyr::select(model, everything()) %>%
    arrange(desc(Accuracy_l1out))
  rownames(res) <- NULL
  return(res)
}
