#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param plim PARAM_DESCRIPTION, Default: 0.01
#' @param opt PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION, Default: 'Condition'
#' @param nfolds PARAM_DESCRIPTION, Default: 0
#' @param xgboost_params PARAM_DESCRIPTION, Default: xgboost_params
#' @param catboost_params PARAM_DESCRIPTION, Default: catboost_params
#' @param randomforest_params PARAM_DESCRIPTION, Default: randomforest_params
#' @param do_smote PARAM_DESCRIPTION, Default: FALSE
#' @param smote_params PARAM_DESCRIPTION, Default: smote_params
#' @param ensemble_param PARAM_DESCRIPTION, Default: 'BalancedAccuracy_l1out'
#' @param ensemble_minval PARAM_DESCRIPTION, Default: 0
#' @param ensemble_1knn PARAM_DESCRIPTION, Default: FALSE
#' @param do_ensemble_probs PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[readr]{write_delim}}
#' @rdname makeAllModels
#' @export 
#' @importFrom readr write_tsv
makeAllModels <- function(datasc, plim=0.01, opt, name="Condition", nfolds=0,
                          xgboost_params = xgboost_params,
                          catboost_params = catboost_params,
                          randomforest_params = randomforest_params,
                          do_smote=FALSE,
                          smote_params=smote_params,
                          ensemble_param = "BalancedAccuracy_l1out",
                          ensemble_minval = 0,
                          ensemble_1knn = FALSE,
                          do_ensemble_probs=TRUE){
  levs <- datasc %>% pull(class) %>% as.factor %>% levels
  # Select features
  if(length(levs)==2){
    compsig <- get_signif_components(datasc, levs)
  }else{
    compsig <- get_signif_components_multiclass(datasc, levs)
  }

  tryCatch({readr::write_tsv(compsig, file=paste0(opt$out, "significant_PCAcomponents_", name,".tsv"))},
           error = function(x){print("ERROR writting sig Components")})

  varnames <- c(compsig$var[compsig$pval <= plim])
  if(length(varnames) < 2){
    varnames <- compsig %>% arrange(pval) %>% head(2) %>% pull(var)
  }

  if(nfolds == 0){
    folds <- c() ## leave 1 out
  }else{
    folds <- createFolds(datasc$class, k = nfolds, list = TRUE, returnTrain = FALSE)
  }


  if(length(levs)==2){
    res_glms <- make_glm_l1o(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params)
  }else{
    res_glms <- make_glm_l1o_multiclass(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params)
  }
  cat("-- GLM finished\n")
  res_svm_lin <- make_svm_l1o(datasc, levs, varnames, kernel="linear", folds = folds, do_smote = do_smote, smote_params = smote_params, balance_classes = TRUE)
  res_svm_rad <- make_svm_l1o(datasc, levs, varnames, kernel="radial", folds = folds, do_smote = do_smote, smote_params = smote_params, balance_classes = FALSE)
  cat("-- SVMs finished\n")
  res_randfor <- make_randomForest_l1o(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params, randomforest_params = randomforest_params)
  cat("-- RandomForest finished\n")
  res_tree <- make_classifTree_l1o(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params, balance_weights = TRUE)
  cat("-- C5.0 Tree finished\n")
  res_naivebayes <- makeNaiveBayes_l1o(datasc, levs, varnames, SEED=SEED, folds = folds, do_smote = do_smote, smote_params = smote_params)
  cat("-- NaiveBayes finished\n")
  res_knn_l1o <- makeKnn_l1o(datasc, levs, varnames, different_ks=seq(3,11, by=2), folds = folds, do_smote = do_smote, smote_params = smote_params)
  #res_knn_no_l1o <- makeKnn(datasc, levs, varnames, different_ks=seq(1,13, by=2))
  cat("-- KNN finished\n")
  res_kmeans_l1o <- makeKmeans_l1o(datasc, levs, varnames, SEED=SEED, folds = folds, do_smote = do_smote, smote_params = smote_params)
  cat("-- K-Means finished\n")
  res_xgboost <- make_xgboost_l1o(datasc, levs, varnames, xgboost_params = xgboost_params, folds = folds, do_smote = do_smote, smote_params = smote_params)
  cat("-- XGBoost finished\n")
  res_catboost <- make_catboost_l1o(datasc, levs, varnames, catboost_params = catboost_params, folds = folds, do_smote = do_smote, smote_params = smote_params)
  cat("-- CatBoost finished\n")

  modlist <- list("logistic_regression" = res_glms,
                  "SVM-linear"=res_svm_lin,
                  "SVM-radial"=res_svm_rad,
                  "RandomForest"=res_randfor,
                  "C5.0 Tree"=res_tree,
                  "NaiveBayes"=res_naivebayes,
                  "XGBoost"=res_xgboost,
                  "CatBoost"=res_catboost,
                  "KMeans"=res_kmeans_l1o)
  if(length(levs)>2)names(modlist)[[1]] <- "Multinom"
  for(k in names(res_knn_l1o)) modlist[[paste0("KNN-", k)]] <- res_knn_l1o[[k]]
  save(modlist, file=paste0(opt$out, "all_models_", name, ".RData"))

  if(length(levs)==2){
    model_res <- getTableFromConfmatrices(modlist)
  }else{
    model_res <- getTableFromConfmatrices_multiclass(modlist)
  }

  modlist$Ensemble <- make_ensemble_votes(datasc, levs, modlist, model_res,
                                          param = ensemble_param,
                                          min_val = ensemble_minval,
                                          only_1_knn = ensemble_1knn)
  cat("-- Ensemble finished\n")
  if(do_ensemble_probs){
    modlist$Ensemble2 <- make_ensemble_probs(datasc, levs, modlist, model_res, param = ensemble_param, min_val = ensemble_minval, only_1_knn = ensemble_1knn)
  }
  if(length(levs)==2){
    model_res <- getTableFromConfmatrices(modlist)
  }else{
    model_res <- getTableFromConfmatrices_multiclass(modlist)
  }

  write_tsv(model_res, file=paste0(opt$out,"summary_all_models", name, ".tsv"))
  return(list(models=modlist, modummary=model_res, component_pvals=compsig, varnames=varnames))
}
