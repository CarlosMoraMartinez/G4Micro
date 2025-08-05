
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param all_pcas PARAM_DESCRIPTION
#' @param PCs PARAM_DESCRIPTION
#' @param modelo_svm PARAM_DESCRIPTION
#' @param vstdf PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param vars2pca PARAM_DESCRIPTION, Default: c("Condition")
#' @param metadata PARAM_DESCRIPTION
#' @param daares PARAM_DESCRIPTION
#' @param topns PARAM_DESCRIPTION, Default: c(5, 10, 20)
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
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}}
#' @rdname callDoAllModelsFromALLPCAsOriginalVars
#' @export 
#' @importFrom dplyr select mutate filter
callDoAllModelsFromALLPCAsOriginalVars <- function(all_pcas, PCs, modelo_svm, vstdf,
                                                   name, vars2pca=c("Condition"), metadata,
                                                   daares, topns = c(5, 10, 20),
                                                   variable_plim=0.01,
                                                   meta_vars = c() ,
                                                   nfolds = 0,
                                                   xgboost_params = xgboost_params,
                                                   catboost_params = catboost_params,
                                                   randomforest_params = randomforest_params,
                                                   do_smote=FALSE,
                                                   smote_params=list(K=5, dup_size="balance")){

  pcts <- summary(all_pcas[[1]]$pca)$importance[2, PCs]
  pcslope <- pcts[1]/pcts[2]
  beta <- drop(t(modelo_svm$coefs) %*% all_pcas$Condition$pca$x[modelo_svm$index,PCs])
  bslope <- -beta[1]/beta[2]

  rotvals <- all_pcas[[1]]$pca$rotation %>% as.data.frame %>% dplyr::select(all_of(PCs)) %>%
    rownames_to_column("taxon") %>%
    dplyr::mutate(
      score1 = abs(all_pcas[[1]]$pca$rotation[, PCs[1]]),
      score2 = abs(all_pcas[[1]]$pca$rotation[, PCs[2]]),
      score3 = abs(all_pcas[[1]]$pca$rotation[, PCs[1]]) + abs(all_pcas[[1]]$pca$rotation[, PCs[2]]),
      score4 = pcslope*abs(all_pcas[[1]]$pca$rotation[, PCs[1]]) + abs(all_pcas[[1]]$pca$rotation[, PCs[2]]),
      score5 = bslope*abs(all_pcas[[1]]$pca$rotation[, PCs[1]]) + abs(all_pcas[[1]]$pca$rotation[, PCs[2]]),
      daa_pval  = -10*log10(daares$padj[match(taxon, daares$taxon)])
    )
  names(rotvals)[!names(rotvals) %in% c("taxon", PCs)] <- c(paste0(PCs[1], " score"),
                                                            paste0(PCs[2], " score"),
                                                            paste0(PCs, collapse="+"),
                                                            paste0(PCs[1], " and ", PCs[2], " combined 2"),
                                                            paste0(PCs[1], " and ", PCs[2], " combined "),
                                                            "DESeq pval"
  )
  modresults <- list()
  for(score in names(rotvals)[!names(rotvals) %in% c("taxon", PCs)]){
    for(topn in topns){
      print(paste0("Fitting models to ", score, '_', topn))
      toptaxa <- rotvals[order(rotvals[, score], decreasing = T), "taxon"][1:topn]
      df2pred <- vstdf %>% dplyr::filter(gene %in% toptaxa) %>% column_to_rownames("gene") %>% as.matrix %>% t %>%
        as.data.frame() %>% rownames_to_column("sample") %>%
        dplyr::mutate(class=unlist(metadata[match(sample, metadata$sampleID), vars2pca[1]]))
      names(df2pred) <- gsub("[\\.\\-\\[\\]()]", "", names(df2pred), perl=T)
      modresults[[paste0(score, ' top ', as.character(topn))]] <- makeAllModels(df2pred, plim=1, opt, name= paste0(name, "_modsIndBacs_", score, "_top", topn),
                                                                                nfolds = nfolds,
                                                                                xgboost_params = xgboost_params,
                                                                                catboost_params = catboost_params,
                                                                                randomforest_params = randomforest_params,
                                                                                do_smote = do_smote, smote_params = smote_params)
      modresults[[paste0(score, ' top ', as.character(topn))]]$taxa <- toptaxa

    }
  }##make models
  print("Merging models")
  modall_table <- map(names(modresults), \(x){
    res <- modresults[[x]]$modummary %>%
      dplyr::mutate(sel_method=x,
                    varsused=paste0(modresults[[x]]$taxa, collapse="|"))

  }) %>% bind_rows()

  return(list(fullresults=modresults, allmodsum=modall_table))
}
