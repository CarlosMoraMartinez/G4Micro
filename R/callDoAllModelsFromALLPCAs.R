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
