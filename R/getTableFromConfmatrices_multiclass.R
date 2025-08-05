#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param modlist PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
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
