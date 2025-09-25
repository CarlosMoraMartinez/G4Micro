#' @title Train and Evaluate Models from PCA Results
#' @description Trains multiple machine learning models using PCA-reduced data and evaluates their performance. Also generates SVM plots for linear and radial kernels.
#' @param all_pcas List of PCA objects as output from prcomp or similar, where the first element must contain `$pca$x`.
#' @param name Character string, prefix for naming the output files.
#' @param metadata Data frame with sample metadata. Must contain column `sampleID` and variables to join and model.
#' @param vars2pca Character vector of metadata variable names to use as class labels for modeling, Default: c("Condition").
#' @param levs2predict Levels of the variable to predict. If empty, gets them automatically in alphabetical order. Default: c()
#' @param variable_plim P-value threshold for variable preselection, based on logistic regressions, Default: 0.01.
#' @param meta_vars Character vector of additional metadata variables to merge into the modeling dataset, Default: c().
#' @param nfolds Integer, number of folds for cross-validation (0 means leave-one-out), Default: 0.
#' @param xgboost_params List of parameters for the XGBoost model, Default: xgboost_params_default.
#' @param catboost_params List of parameters for the CatBoost model, Default: catboost_params_default.
#' @param randomforest_params List of parameters for the Random Forest model, Default: randomforest_params_default.
#' @param do_smote Logical, whether to apply SMOTE oversampling, Default: FALSE.
#' @param smote_params List of parameters for SMOTE oversampling, Default: list(K = 5, dup_size = "balance").
#' @param opt A list of options, must contain an output path in `opt$out`.
#' @return A list containing model summaries and generated plots, including trained models and SVM visualizations.
#' @details This function prepares a dataset from PCA results and associated metadata, selects variables based on logistic regression (for binomial data), trains various machine learning models (SVM, XGBoost, CatBoost, Random Forest, KNN, K-Means, C5.0 Tree, NaiveBayes, GLM), and optionally applies SMOTE. It returns a summary and generates plots for linear and radial SVMs.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  results <- callDoAllModelsFromALLPCAs(all_pcas, name = "example", metadata = meta_df)
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}},
#'  \code{\link[dplyr]{filter}},
#'  \code{\link[dplyr]{select}},
#'  \code{\link[dplyr]{inner_join}},
#'  \code{\link[tibble]{rownames_to_column}},
#'  \code{\link[dplyr]{all_of}},
#'  \code{\link[dplyr]{join_by}},
#'  \code{\link{makeAllModels}},
#'  \code{\link{plotSVM}}
#' @rdname callDoAllModelsFromALLPCAs
#' @export
#' @importFrom dplyr mutate filter select
callDoAllModelsFromALLPCAs <- function(all_pcas, name, metadata, vars2pca=c("Condition"),
                                       levs2predict = c(),
                                       variable_plim=0.01,
                                       meta_vars = c(),
                                       nfolds = 0,
                                       xgboost_params = xgboost_params_default,
                                       catboost_params = catboost_params_default,
                                       randomforest_params = randomforest_params_default,
                                       do_smote=FALSE,
                                       smote_params=smote_params_default, opt){
  datasc <- all_pcas[[1]]$pca$x %>%
    as.data.frame %>%
    rownames_to_column("sample") %>%
    dplyr::mutate(class=unlist(metadata[match(sample, metadata$sampleID), vars2pca[1]])) %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::mutate(class=factor(class))
  if(length(meta_vars) > 0){
    meta_filt <- metadata %>% dplyr::select(sampleID, all_of(meta_vars))
    byy <- join_by(sample == sampleID)
    datasc <- datasc %>% dplyr::inner_join(meta_filt, by=byy)

  }
  allmodssumm <- makeAllModels(datasc, plim=variable_plim, opt,
                               name= name, nfolds = nfolds,
                               levs2predict = levs2predict,
                               xgboost_params = xgboost_params,
                               catboost_params = catboost_params,
                               randomforest_params = randomforest_params,
                               do_smote = do_smote, smote_params = smote_params,
                               do_ensemble_probs = FALSE)

  modelo_svm <- allmodssumm$models$`SVM-linear`$mod_noscale
  allmodssumm$plot_svm_rad <-plotSVM(modelo_svm, datasc, allmodssumm$varnames,
                                     opt, paste0(name, "_linear"))

  modelo_svm <- allmodssumm$models$`SVM-radial`$mod_noscale
  allmodssumm$plot_svm_rad <- plotSVM(modelo_svm, datasc, allmodssumm$varnames,
                                      opt, paste0(name, "_radial"))

  return(allmodssumm)
}
