#' @title Build and Evaluate Multiple Machine Learning Models with K-Fold Cross-Validation or Leave-One-Out Strategy
#' @description This function builds a set of machine learning models using leave-one-out or k-fold cross-validation.
#' It selects significant PCA components, applies various classifiers, and optionally ensembles them.
#' @param datasc A data frame with a column named `class` containing the labels.
#' @param plim p-value threshold to select significant variables. Default: 0.01
#' @param opt A list of options, must contain an output path in `opt$out`.
#' @param name Name for the output files and modeling group. Default: 'Condition'
#' @param nfolds Number of folds for cross-validation. If 0, leave-one-out is used. Default: 0
#' @param levs2predict Levels of the variable to predict. If empty, gets them automatically in alphabetical order. Default: c()
#' @param xgboost_params Parameters for the xgboost model. Default: xgboost_params_default
#' @param catboost_params Parameters for the catboost model. Default: catboost_params_default
#' @param randomforest_params Parameters for the random forest model. Default: randomforest_params_default
#' @param do_smote Whether to apply SMOTE to balance classes. Default: FALSE
#' @param smote_params Parameters for SMOTE. Default: smote_params_default
#' @param ensemble_param Performance metric used for selecting models in the ensemble. Default: 'BalancedAccuracy_l1out'
#' @param ensemble_minval Minimum performance required for inclusion in ensemble. Default: 0
#' @param ensemble_1knn Whether to include only 1-NN in the ensemble. Default: FALSE
#' @param do_ensemble_probs Whether to build a second ensemble using probability-based voting. Default: TRUE
#' @return A list containing: \code{models}, \code{modummary} (performance summary), \code{component_pvals} (significant PCA components), and \code{varnames} (selected variable names).
#' @details This function automates model training and evaluation using various classifiers (GLM, SVM, Random Forest, etc.) and ensembling.
#' It is designed for multi-class or binary classification.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   result <- makeAllModels(datasc, opt = list(out = "output_dir/"))
#' }
#' }
#' @seealso
#'  \code{\link[readr]{write_tsv}},
#'  \code{\link{get_signif_components}},
#'  \code{\link{get_signif_components_multiclass}},
#'  \code{\link{createFolds}},
#'  \code{\link{make_glm_l1o}},
#'  \code{\link{make_glm_l1o_multiclass}},
#'  \code{\link{make_svm_l1o}},
#'  \code{\link{make_randomForest_l1o}},
#'  \code{\link{make_classifTree_l1o}},
#'  \code{\link{makeNaiveBayes_l1o}},
#'  \code{\link{makeKnn_l1o}},
#'  \code{\link{makeKmeans_l1o}},
#'  \code{\link{make_xgboost_l1o}},
#'  \code{\link{make_catboost_l1o}},
#'  \code{\link{getTableFromConfmatrices}},
#'  \code{\link{getTableFromConfmatrices_multiclass}},
#'  \code{\link{make_ensemble_votes}},
#'  \code{\link{make_ensemble_probs}}
#' @rdname makeAllModels
#' @export
#' @importFrom readr write_tsv
#' @importFrom caret createFolds
makeAllModels <- function(datasc, plim=0.01, opt, name="Condition", nfolds=0,
                          levs2predict = c(),
                          xgboost_params = xgboost_params_default,
                          catboost_params = catboost_params_default,
                          randomforest_params = randomforest_params_default,
                          do_smote=FALSE,
                          smote_params=smote_params_default,
                          ensemble_param = "BalancedAccuracy_l1out",
                          ensemble_minval = 0,
                          ensemble_1knn = FALSE,
                          do_ensemble_probs=TRUE){
  if(length(levs2predict) == 0){
    levs <- datasc %>% pull(class) %>% as.factor %>% levels
  }else{
    levs <- levs2predict
  }
  datasc <- datasc %>% dplyr::mutate(class=factor(class, levels=levs))
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
    folds <- caret::createFolds(datasc$class, k = nfolds, list = TRUE, returnTrain = FALSE)
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
  res_tree <- make_classifTree_l1o(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params, balance_weights = FALSE)
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
