#' @title Naive Bayes Classifier with General Cross-Validation and Optional SMOTE Oversampling
#' @description
#' Trains and evaluates a Naive Bayes classifier using general cross-validation specified by \code{folds}.
#' Optionally applies SMOTE oversampling to training folds to address class imbalance.
#' Returns confusion matrices, predictions, predicted probabilities, and ROC AUC metrics.
#'
#' @param datasc A data frame containing the dataset, including a factor column named \code{class} for the target variable,
#' and a \code{sample} column. Predictor variables are selected by \code{varnames}.
#' @param levs A character vector specifying the factor levels of the target variable \code{class}. The second element is used as the positive class.
#' @param varnames A character vector of predictor variable names to use in the model.
#' @param SEED Integer seed for random number generation to ensure reproducibility. Default is 123.
#' @param folds Integer vector or list specifying the indices of samples used as test sets in cross-validation.
#' If empty, performs leave-one-out cross-validation (each sample left out once).
#' @param do_smote Logical indicating whether to apply SMOTE oversampling on training data. Default is \code{FALSE}.
#' @param smote_params Named list with SMOTE parameters:
#' \itemize{
#'   \item \code{K}: Number of nearest neighbors (default 5).
#'   \item \code{dup_size}: Oversampling amount, can be \code{"balance"} or numeric (default \code{"balance"}).
#' }
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{\code{confmat}}{Confusion matrix from cross-validated predictions.}
#'   \item{\code{confmat_no_l1o}}{Confusion matrix from the model trained on the entire dataset without cross-validation.}
#'   \item{\code{preds}}{Vector of predicted classes from cross-validation.}
#'   \item{\code{pred_probs}}{Predicted probabilities for the positive class from cross-validation.}
#'   \item{\code{preds_no_l1o}}{Predicted classes from the model trained on the entire dataset.}
#'   \item{\code{mod}}{Naive Bayes model trained on the entire dataset.}
#'   \item{\code{roc_obj_no_l1o}}{Currently \code{NULL}, placeholder for ROC object from non-cross-validated model.}
#'   \item{\code{roc_auc_no_l1o}}{Currently \code{NULL}, placeholder for ROC AUC value from non-cross-validated model.}
#'   \item{\code{roc_obj}}{ROC object from cross-validation predictions (binary or multiclass).}
#'   \item{\code{roc_auc}}{Numeric ROC AUC value from cross-validation.}
#' }
#'
#' @details
#' This function performs general cross-validation for a Naive Bayes classifier using the \code{e1071} package.
#' If \code{folds} is empty, it defaults to leave-one-out cross-validation.
#' When \code{do_smote} is \code{TRUE}, SMOTE oversampling is applied on the training fold to improve class balance, using \code{SmoteClassif}.
#' ROC and AUC metrics are computed differently for binary and multiclass classification.
#' The function also fits a Naive Bayes model on the entire dataset for comparison.
#'
#' The user must ensure that the \code{datasc} data frame contains the required columns and that \code{levs} corresponds to levels in \code{datasc$class}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'
#'   # Example: use iris dataset (binary subset)
#'   data(iris)
#'   iris_sub <- iris %>% filter(Species != "setosa")
#'   iris_sub$class <- factor(iris_sub$Species)
#'   iris_sub$sample <- 1:nrow(iris_sub)
#'   levs <- levels(iris_sub$class)
#'   varnames <- colnames(iris_sub)[1:4]
#'
#'   results <- makeNaiveBayes_l1o(
#'     datasc = iris_sub,
#'     levs = levs,
#'     varnames = varnames,
#'     SEED = 42,
#'     folds = c(),         # leave-one-out
#'     do_smote = FALSE
#'   )
#'
#'   print(results$confmat)
#'   print(results$roc_auc)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{select}},
#' \code{\link[e1071]{naiveBayes}},
#' \code{\link[caret]{confusionMatrix}},
#' \code{\link[pROC]{roc}},
#' \code{\link[UBL]{SmoteClassif}},
#' \code{\link[purrr]{map}},
#' \code{\link[rlang]{sym}}
#'
#' @rdname makeNaiveBayes_l1o
#' @export
#' @importFrom dplyr select
#' @importFrom e1071 naiveBayes
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc multiclass.roc
#' @importFrom UBL SmoteClassif
#' @importFrom purrr map
#' @importFrom rlang sym
makeNaiveBayes_l1o <- function(datasc, levs, varnames,
                               SEED=123, folds=c(),
                               do_smote=FALSE,
                               smote_params=list(K=5, dup_size="balance")){
  set.seed(SEED)

  predict_bayes1 <- factor()
  predict_bayes1_probs <- list()
  df <- datasc %>% dplyr::select(-class, -sample) %>% dplyr::select(all_of(varnames))
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]

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
    mod_bayes2 <- naiveBayes(train_df, train_labels, laplace = 0)
    predict_bayes1 <- c(predict_bayes1, predict(mod_bayes2, test_df))
    predict_bayes1_probs[[i]] <- predict(mod_bayes2, test_df, type="raw")

  }
  confusionMatrix_bayes1 <- confusionMatrix(predict_bayes1, datasc$class,
                                            positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict_bayes1_probs, \(xx) xx %>% as.data.frame) %>%
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- map(predict_bayes1_probs, \(xx) xx %>% as.data.frame) %>%
      bind_rows
    #probs_vector2 <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>%
    #   map2(datasc$class, \(x, nn) x[as.character(nn)]) %>% unlist
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
  }
  modwithall <- e1071::naiveBayes(df, datasc$class, laplace = 0)
  predict_bayes2 <- predict(modwithall, df)
  confMatrix_bayes2_nol1o <- confusionMatrix(predict_bayes2, datasc$class,
                                             positive = levs[2])
  return(list(confmat=confusionMatrix_bayes1,
              confmat_no_l1o=confMatrix_bayes2_nol1o,
              preds=predict_bayes1,
              pred_probs=probs_vector,
              preds_no_l1o=predict_bayes2,
              mod=modwithall,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc))
  )
}
