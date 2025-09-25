#' @title Train and Evaluate an SVM Classifier with Cross-Validation
#' @description
#' This function trains a Support Vector Machine (SVM) classifier using the `e1071` package and evaluates it using cross-validation.
#' It supports both leave-one-out and custom-fold validation, optional SMOTE oversampling, and returns detailed performance metrics
#' including confusion matrices and ROC AUC.
#'
#' @param datasc A data frame containing the input data. Must include a `class` column (as factor) and the specified feature columns.
#' @param levs A character vector specifying the class levels (e.g., `c("control", "case")`).
#' @param varnames A character vector with the names of the variables (features) to be used for training.
#' @param kernel The kernel to be used by the SVM model. Default: `'linear'`.
#' @param SEED An integer used to set the random seed for reproducibility. Default: `123`.
#' @param folds A vector of indices indicating which samples to leave out in each iteration. If empty, leave-one-out is used. Default: `c()`.
#' @param do_smote Logical, whether to apply SMOTE to balance the training data. Default: `FALSE`.
#' @param smote_params A list of parameters passed to `SmoteClassif`, including `K` (number of neighbors) and `dup_size` (oversampling factor). Default: `list(K = 5, dup_size = 2)`.
#' @param balance_classes Logical, whether to balance class weights in the SVM. Default: `TRUE`.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item `confmat`: Confusion matrix from cross-validation predictions.
#'   \item `confmat_no_l1o`: Confusion matrix from model trained on all data.
#'   \item `mod`: Final SVM model trained on all data.
#'   \item `preds`: Predictions from cross-validation.
#'   \item `pred_probs`: Probabilities from cross-validation predictions.
#'   \item `pred_probs_obj`: Raw probability prediction objects from cross-validation.
#'   \item `preds_no_l1o`: Predictions from model trained on all data.
#'   \item `mod_noscale`: SVM model trained on first 2 features without scaling.
#'   \item `preds_noscale`: Predictions from `mod_noscale`.
#'   \item `confmat_noscale`: Confusion matrix for `preds_noscale`.
#'   \item `roc_obj`: ROC object from cross-validation.
#'   \item `roc_auc`: AUC value from `roc_obj`.
#'   \item `roc_obj_no_l1o`: Placeholder, always NULL.
#'   \item `roc_auc_no_l1o`: Placeholder, always NULL.
#' }
#'
#' @details
#' The function allows for either leave-one-out or custom-fold cross-validation using SVM classifiers.
#'  It optionally applies SMOTE oversampling and computes ROC AUC using either binary or multiclass evaluation.
#'  SVMs are trained using the `e1071` package with optional class weighting.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  result <- make_svm_l1o(datasc = my_data,
#'                         levs = c("control", "case"),
#'                         varnames = c("gene1", "gene2", "gene3"),
#'                         kernel = "linear")
#'  print(result$confmat)
#'  plot(result$roc_obj)
#'  }
#' }
#'
#' @seealso
#'  \code{\link[dplyr]{select}},
#'  \code{\link[dplyr]{all_of}},
#'  \code{\link[e1071]{svm}},
#'  \code{\link[UBL]{SmoteClassif}},
#'  \code{\link[pROC]{roc}},
#'  \code{\link[pROC]{multiclass.roc}},
#'  \code{\link[caret]{confusionMatrix}},
#'  \code{\link[purrr]{map}},
#'  \code{\link[dplyr]{bind_rows}},
#'  \code{\link[rlang]{sym}},
#'  \code{\link[stats]{as.formula}},
#'  \code{\link[base]{set.seed}},
#'  \code{\link[base]{factor}},
#'  \code{\link[base]{rep}},
#'  \code{\link[base]{paste}},
#'  \code{\link[base]{paste0}}
#'
#' @rdname make_svm_l1o
#' @export
#' @importFrom dplyr select all_of bind_rows
#' @importFrom e1071 svm
#' @importFrom UBL SmoteClassif
#' @importFrom pROC roc multiclass.roc
#' @importFrom caret confusionMatrix
#' @importFrom purrr map
#' @importFrom rlang sym
#' @importFrom stats as.formula
make_svm_l1o <- function(datasc, levs, varnames, kernel="linear", SEED=123, folds=c(),
                         do_smote=FALSE,
                         smote_params=smote_params_default,
                         balance_classes=TRUE){
  datasc$class <- factor(datasc$class)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  predict1 <- factor(levels=levs)
  predict1_probs <- list()

  if(balance_classes & ! do_smote){
    class_weight <- "inverse"
  }else{
    class_weight <- rep(1, length(levs))
    names(class_weight) <- levs
  }

  if(length(folds)==0){
    folds <- 1:nrow(datasc)
    reorder_samples <- FALSE
  }else{
    reorder_samples <- TRUE
  }
  set.seed(SEED)
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
      train_df <- smoteData %>% dplyr::select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
    }else{
      smoteData = NULL
    }

    mod <- e1071::svm(x = train_df, y = train_labels, scale=TRUE, kernel=kernel,
                      class.weights = class_weight,
                      probability = TRUE)
    predict1 <- c(predict1, predict(mod, test_df))
    if(length(i) > 1){
      for(ii in i)predict1_probs[[ii]] <- predict(mod, test_df[as.character(ii), ], probability = TRUE)
    }else{
      predict1_probs[[i]] <- predict(mod, test_df, probability = TRUE)
    }

  }
  if(reorder_samples){
    predict1 <- predict1[as.character(1:nrow(datasc))]
  }
  confmat1 <- confusionMatrix(predict1, datasc$class, positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>%
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- pROC::roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>%
      bind_rows #%>% pull(!!sym(levs[2]))
    roc1 <- pROC::multiclass.roc(response=datasc$class, predictor=probs_vector)
    roc_auc <- as.numeric(roc1$auc)
  }

  mod_all <- e1071::svm(x = df, y = datasc$class, scale=TRUE, kernel=kernel,
                        probability = TRUE, class.weights = class_weight)
  predict2 <- predict(mod_all, df)
  predict2_probs <- predict(mod_all, df, probability = TRUE)
  confmat2 <- confusionMatrix(predict2, datasc$class, positive = levs[2])

  mod_all_noscale <- e1071::svm(x = df[, varnames[1:2]], y = datasc$class, scale=FALSE,
                                kernel=kernel, class.weights = "inverse")
  predict2_noscale <- predict(mod_all_noscale, df[, varnames[1:2]])
  confmat2_noscale <- confusionMatrix(predict2_noscale, datasc$class, positive = levs[2])

  return(list(confmat=confmat1,
              confmat_no_l1o=confmat2,
              mod=mod_all,
              preds=predict1,
              pred_probs = probs_vector,
              pred_probs_obj = predict1_probs,
              preds_no_l1o=predict2,
              mod_noscale=mod_all_noscale,
              preds_noscale=predict2_noscale,
              confmat_noscale=confmat2_noscale,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc))
  )
}
