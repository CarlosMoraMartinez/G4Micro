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
#' @rdname getTableFromConfmatrices
#' @export 
#' @importFrom dplyr mutate select
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
