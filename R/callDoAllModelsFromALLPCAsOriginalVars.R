
#' @title Train models using selected PCA rotations and DA results
#' @description
#' This function selects important taxa from PCA loadings and differential abundance results (e.g., DESeq2),
#' and trains predictive models (XGBoost, CatBoost, Random Forest, optionally with SMOTE) using subsets
#' of features. It supports general K-fold cross-validation.
#'
#' The input data must contain a column named `'sample'` to match sample-level metadata.
#'
#' @param all_pcas A list of PCA results, each with an entry `$pca` (result from `prcomp`).
#' @param PCs Vector with the two principal components to use for ranking taxa (e.g., c(1,2)).
#' @param modelo_svm SVM model trained on PCA-transformed data; used to compute class separation direction.
#' @param vstdf Data frame with variance-stabilized abundance data. Must contain a column named `gene`.
#' @param name Character name used as prefix for saved model names.
#' @param vars2pca Character vector of variable names from metadata used to train models. Default: `"Condition"`.
#' @param metadata Data frame with sample-level metadata. Must contain a `sampleID` column.
#' @param daares Data frame with differential abundance results. Must contain columns `taxon` and `padj`.
#' @param topns Numeric vector of how many top-ranked taxa to use for model training. Default: c(5, 10, 20).
#' @param variable_plim Numeric threshold for filtering variables by p-value. Default: 0.01.
#' @param meta_vars Character vector of additional metadata variables to include in modeling. Default: empty.
#' @param nfolds Number of folds to use for cross-validation. Use `nfolds = 0` for no cross-validation. Default: 0.
#' @param xgboost_params List of parameters for XGBoost training. Default: `xgboost_params_default`.
#'   A typical default might be:
#'   ```
#'   list(
#'     learning_rate = 0.3, max_depth = 2, nrounds = 30, min_child_weight = 1,
#'     subsample = 1, colsample_bytree = 0.6, reg_lambda = 1, reg_alpha = 0,
#'     nthread = 1, objective = "multi:softprob", num_class = 4, balance_weights = TRUE
#'   )
#'   ```
#' @param catboost_params List of parameters for CatBoost training. Default: `catboost_params_default`.
#' @param randomforest_params List of parameters for Random Forest training. Default: `randomforest_params_default`.
#' @param do_smote Logical, whether to apply SMOTE for class balancing before training. Default: FALSE.
#' @param smote_params List with SMOTE configuration. Default: `smote_params_default`.
#'
#' @return A list with:
#' \itemize{
#'   \item `fullresults`: A named list of model objects for each (score, topn) combination.
#'   \item `allmodsum`: A summary data frame of model performances and used taxa.
#' }
#'
#' @details
#' Taxa are ranked based on their contribution to PCA components and optionally the direction of the SVM decision boundary.
#' Models are trained on top-N taxa using different scoring strategies.
#'
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
                                                   xgboost_params = xgboost_params_default,
                                                   catboost_params = catboost_params_default,
                                                   randomforest_params = randomforest_params_default,
                                                   do_smote=FALSE,
                                                   smote_params=smote_params_default){

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
