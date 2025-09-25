#' @title Classification Tree with General Cross-Validation and Optional SMOTE Oversampling
#' @description
#' Trains and evaluates a classification tree using C5.0 with general cross-validation defined by \code{folds}.
#' Optionally applies SMOTE oversampling to training data folds to handle class imbalance.
#' Returns confusion matrices, predictions, predicted probabilities, and ROC AUC metrics.
#'
#' @param datasc A data frame containing the dataset. Must include a factor column \code{class} with the target labels,
#' and a \code{sample} column. Predictor variables are selected by \code{varnames}.
#' @param levs A character vector specifying the factor levels of the \code{class} variable. The second element is used as the positive class.
#' @param varnames A character vector of predictor variable names to use in the model.
#' @param folds Integer vector specifying indices of samples to be used as test sets in cross-validation. If empty, leave-one-out cross-validation is performed.
#' @param balance_weights Logical indicating whether to apply class weights to balance classes. Default: \code{TRUE}.
#' Note: currently not implemented for C5.0 (warning issued).
#' @param do_smote Logical indicating whether to apply SMOTE oversampling on training folds. Default: \code{FALSE}.
#' @param smote_params Named list with SMOTE parameters: \code{K} (nearest neighbors, default 5), \code{dup_size} (oversampling amount, default \code{"balance"}).
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{confmat}}{Confusion matrix from cross-validated predictions.}
#'   \item{\code{mod}}{C5.0 classification tree model trained on the full dataset.}
#'   \item{\code{preds}}{Predicted classes from cross-validation.}
#'   \item{\code{pred_probs}}{Predicted probabilities for the positive class from cross-validation.}
#'   \item{\code{preds_no_l1o}}{Predicted classes from model trained on the full dataset.}
#'   \item{\code{confmat_no_l1o}}{Confusion matrix from the full dataset model predictions.}
#'   \item{\code{roc_obj_no_l1o}}{Currently \code{NULL}, placeholder for ROC object from the full dataset model.}
#'   \item{\code{roc_auc_no_l1o}}{Currently \code{NULL}, placeholder for ROC AUC value from the full dataset model.}
#'   \item{\code{roc_obj}}{ROC curve object from cross-validation predictions (binary or multiclass).}
#'   \item{\code{roc_auc}}{Numeric ROC AUC value from cross-validation predictions.}
#' }
#'
#' @details
#' This function performs classification using C5.0 trees with cross-validation defined by \code{folds}.
#' If \code{folds} is empty, leave-one-out cross-validation is performed.
#' When \code{do_smote} is \code{TRUE}, SMOTE oversampling is applied on each training fold using parameters defined in \code{smote_params}.
#' Class weights balancing is currently not implemented due to C5.0 limitations.
#' ROC/AUC is computed differently depending on whether the classification is binary or multiclass.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'
#'   data(iris)
#'   iris_sub <- iris %>% filter(Species != "setosa")
#'   iris_sub$class <- factor(iris_sub$Species)
#'   iris_sub$sample <- 1:nrow(iris_sub)
#'   levs <- levels(iris_sub$class)
#'   varnames <- colnames(iris_sub)[1:4]
#'
#'   res <- make_classifTree_l1o(
#'     datasc = iris_sub,
#'     levs = levs,
#'     varnames = varnames,
#'     folds = c(),           # leave-one-out CV
#'     balance_weights = FALSE,
#'     do_smote = FALSE
#'   )
#'   print(res$confmat)
#'   print(res$roc_auc)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{select}},
#'  \code{\link[C50]{C5.0}},
#'  \code{\link[caret]{confusionMatrix}},
#'  \code{\link[pROC]{roc}},
#'  \code{\link[UBL]{SmoteClassif}},
#'  \code{\link[purrr]{map}},
#'  \code{\link[rlang]{sym}}
#'
#' @rdname make_classifTree_l1o
#' @export
#' @importFrom dplyr select
#' @importFrom C50 C5.0
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc multiclass.roc
#' @importFrom UBL SmoteClassif
#' @importFrom purrr map
#' @importFrom rlang sym
make_classifTree_l1o <- function(datasc, levs, varnames,
                                 folds=c(),
                                 balance_weights = TRUE,
                                 do_smote=FALSE,
                                 smote_params=list(K=5, dup_size="balance")){
  predict_tree1 <- list()
  predict_probs <- list()
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  if(balance_weights & !do_smote){
    warning("Using weights in C5.0 not implemented")
    classweights <- table(datasc$class)
    sweights <- max(classweights)/classweights[datasc$class]
  }else{
    sweights <- rep(1, nrow(datasc))
  }
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
    reorder_samples <- FALSE
  }else{
    reorder_samples <- TRUE
  }
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    #train_weighs <- sweights[-i]

    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]

    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df,
                                C.perc = smote_params$dup_size,
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
    }else{
      smoteData = NULL
    }

    mod_tree1 <- C5.0(train_df, train_labels, trials = 20) # , weights=train_weighs # makes it crash!

    if(length(i) > 1){
      for(ii in i) predict_tree1[[ii]] <- predict(mod_tree1, test_df[as.character(ii), ])
      for(ii in i) predict_probs[[ii]] <- predict(mod_tree1, test_df[as.character(ii), ], type = "prob")
    }else{
      predict_tree1 <- c(predict_tree1, predict(mod_tree1, test_df))
      predict_probs[[i]] <- predict(mod_tree1, test_df, type = "prob")
    }
  }
  if(reorder_samples){
    predict_tree1 <- unlist(predict_tree1)
  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict_probs, \(xx) xx %>% as.data.frame) %>%
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- map(predict_probs, \(xx) xx %>% as.data.frame) %>%
      bind_rows #%>% pull(!!sym(levs[2]))
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
  }
  mod_all <- C5.0(df, datasc$class, trials = 20) # , weights=sweighs
  predict_tree2 <- predict(mod_all, df)
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  return(list(confmat=confmat_tree1,
              mod=mod_all,
              preds=predict_tree1,
              pred_probs=probs_vector,
              preds_no_l1o=predict_tree2,
              confmat_no_l1o=confmat_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}
