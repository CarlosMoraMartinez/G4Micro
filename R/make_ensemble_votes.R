#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param levs PARAM_DESCRIPTION
#' @param modlist PARAM_DESCRIPTION
#' @param model_res PARAM_DESCRIPTION
#' @param param PARAM_DESCRIPTION, Default: 'Kappa_l1out'
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
#'  \code{\link[purrr]{map}}
#' @rdname make_ensemble_votes
#' @export 
#' @importFrom dplyr arrange filter
#' @importFrom purrr map
make_ensemble_votes <- function(datasc, levs, modlist, model_res, param="Kappa_l1out",
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
  preddf <- map(m2use, \(x) tibble( !!x := modlist[[x]]$preds))  %>% bind_cols()
  #names(preddf) <- m2use
  if(prop){
    ponderfac <- model_res[match(m2use, model_res$model), param]
    ponderfac <- (ponderfac - min(ponderfac))/(max(ponderfac) - min(ponderfac)) + 0.1
  }else{
    ponderfac <- rep(1, length(m2use))
  }
  preds <- c()
  votes <- list()
  for(i in 1:nrow(preddf)){
    classcore <- map_vec(levs, \(ll) sum(ponderfac[preddf[i, ] == ll]))
    names(classcore) <- levs
    l1 <- length(preds)
    preds <- c(preds, levs[which.max(classcore)] )
    l2 <- length(preds)
    cat(i, ": L1=", l1, ", L2=", l2, ifelse(l1==l2, " --WARNING--", ""),  "\n")
    votes[[i]] <- classcore
  }
  preds <- factor(preds, levels=levs)
  confmat1 <- confusionMatrix(preds, datasc$class, positive = levs[2])

  if(length(levs) == 2){
    probs <- map_vec(votes, \(x)x[levs[2]]/sum(x) )
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs)
  }else{
    probs <- purrr::map(votes, .f = \(x) x/sum(x)) %>%
      bind_rows %>% as.matrix #%>% pull(!!sym(levs[2]))
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs)
  }
  return(list(confmat=confmat1,
              confmat_no_l1o=NULL,
              mod=NULL,
              preds=preds,
              pred_df=preddf,
              preds_no_l1o=NULL,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}
