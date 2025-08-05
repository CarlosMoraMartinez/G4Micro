#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param levs PARAM_DESCRIPTION
#' @param modlist PARAM_DESCRIPTION
#' @param model_res PARAM_DESCRIPTION
#' @param param PARAM_DESCRIPTION, Default: 'BalancedAccuracy_l1out'
#' @param min_val PARAM_DESCRIPTION, Default: 0
#' @param prop PARAM_DESCRIPTION, Default: TRUE
#' @param only_1_knn PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{filter}}
#' @rdname make_ensemble_probs
#' @export 
#' @importFrom dplyr arrange filter
make_ensemble_probs <- function(datasc, levs, modlist, model_res, param="BalancedAccuracy_l1out", #"Kappa_l1out",
                                min_val=0, prop=TRUE,
                                only_1_knn=FALSE){ # 0.65
  remove_knn <- model_res %>%
    dplyr::arrange(desc(!!sym(param))) %>%
    dplyr::filter(grepl("KNN", model)) %>% pull(model)
  remove_knn <- remove_knn[2:length(remove_knn)]
  m2use <- model_res %>%
    dplyr::filter(!!sym(param) >= min_val) %>%
    dplyr::filter(! (model %in% remove_knn & only_1_knn)) %>%
    pull(model)
  preddf <- map(m2use, \(x) modlist[[x]]$pred_probs)  %>% bind_cols()
  names(preddf) <- m2use
  if(prop){
    ponderfac <- model_res[match(m2use, model_res$model), param]
    ponderfac <- (ponderfac - min(ponderfac))/(max(ponderfac) - min(ponderfac)) + 0.1

  }else{
    ponderfac <- rep(1, length(m2use))
  }
  preds <- c()
  avg_probs <- c()
  for(i in 1:nrow(preddf)){
    classcore <-sum(preddf[i, ]*ponderfac)/sum(ponderfac)
    avg_probs <- c(avg_probs, classcore)
    preds <- c(preds, levs[as.integer(round(classcore))+1] )
  }
  preds <- factor(preds, levels=levs)
  confmat1 <- confusionMatrix(preds, datasc$class, positive = levs[2])

  roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=avg_probs)

  return(list(confmat=confmat1,
              confmat_no_l1o=NULL,
              mod=NULL,
              preds=preds,
              pred_probs=avg_probs,
              pred_df=preddf,
              preds_no_l1o=NULL,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}
