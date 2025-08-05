#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param all_pcas PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param metadata PARAM_DESCRIPTION
#' @param vars2pca PARAM_DESCRIPTION, Default: c("Condition")
#' @param variable_plim PARAM_DESCRIPTION, Default: 0.01
#' @param meta_vars PARAM_DESCRIPTION, Default: c()
#' @param nfolds PARAM_DESCRIPTION, Default: 0
#' @param xgboost_params PARAM_DESCRIPTION, Default: xgboost_params
#' @param catboost_params PARAM_DESCRIPTION, Default: catboost_params
#' @param randomforest_params PARAM_DESCRIPTION, Default: randomforest_params
#' @param do_smote PARAM_DESCRIPTION, Default: FALSE
#' @param smote_params PARAM_DESCRIPTION, Default: list(K = 5, dup_size = "balance")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}}, \code{\link[dplyr]{select}}
#' @rdname callDoAllModelsFromALLPCAs
#' @export 
#' @importFrom dplyr mutate filter select
callDoAllModelsFromALLPCAs <- function(all_pcas, name, metadata, vars2pca=c("Condition"),
                                       variable_plim=0.01,
                                       meta_vars = c(),
                                       nfolds = 0,
                                       xgboost_params = xgboost_params,
                                       catboost_params = catboost_params,
                                       randomforest_params = randomforest_params,
                                       do_smote=FALSE,
                                       smote_params=list(K=5, dup_size="balance")){
  datasc <- all_pcas[[1]]$pca$x %>%
    as.data.frame %>%
    rownames_to_column("sample") %>%
    dplyr::mutate(class=unlist(metadata[match(sample, metadata$sampleID), vars2pca[1]])) %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::mutate(class=factor(class))
  if(length(meta_vars) > 0){
    meta_filt <- metadata %>% dplyr::select(sampleID, all_of(meta_vars))
    byy <- join_by(sample == sampleID)
    datasc <- datasc %>% inner_join(meta_filt, by=byy)

  }
  allmodssumm <- makeAllModels(datasc, plim=variable_plim, opt, name= name, nfolds = nfolds,
                               xgboost_params = xgboost_params,
                               catboost_params = catboost_params,
                               randomforest_params = randomforest_params,
                               do_smote = do_smote, smote_params = smote_params,
                               do_ensemble_probs = FALSE)

  modelo_svm <- allmodssumm$models$`SVM-linear`$mod_noscale
  allmodssumm$plot_svm_rad <-plotSVM(modelo_svm, datasc, allmodssumm$varnames, opt, paste0(name, "_linear"))

  modelo_svm <- allmodssumm$models$`SVM-radial`$mod_noscale
  allmodssumm$plot_svm_rad <- plotSVM(modelo_svm, datasc, allmodssumm$varnames, opt, paste0(name, "_radial"))

  return(allmodssumm)
}
