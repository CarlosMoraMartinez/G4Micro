#' @title Multinomial GLM with Custom Fold Cross-Validation and Optional SMOTE
#' @description Trains a multinomial logistic regression model using user-defined folds (e.g., leave-one-out or k-fold) with optional SMOTE resampling to handle class imbalance.
#' @param datasc A data frame containing the classification variable `class`, a `sample` column, and predictor variables.
#' @param levs Character vector with the class levels to use in the factor definition of the response variable.
#' @param varnames Character vector indicating the names of the predictor variables to include in the model.
#' @param folds Integer vector of row indices defining the test folds. If empty (default), performs leave-one-out cross-validation.
#' @param do_smote Logical indicating whether to apply SMOTE resampling to balance the training data, Default: FALSE
#' @param smote_params A list with parameters for SMOTE: `K` (number of neighbors) and `dup_size` (duplication factor), Default: list(K = 5, dup_size = 2)
#' @return A list containing:
#' \itemize{
#'   \item \code{confmat}: Confusion matrix from fold-based prediction.
#'   \item \code{confmat_no_l1o}: Confusion matrix from the full model.
#'   \item \code{mod}: Fitted model using all data.
#'   \item \code{preds}: Class predictions from cross-validation.
#'   \item \code{preds_no_l1o}: Class predictions from the full model.
#'   \item \code{roc_obj}: Multiclass ROC object from fold-based prediction.
#'   \item \code{roc_auc}: AUC from the cross-validation model.
#'   \item \code{roc_obj_no_l1o}: List of ROC objects per class from the full model.
#'   \item \code{roc_auc_no_l1o}: Mean AUC of ROC objects from the full model.
#' }
#' @details
#' This function supports multi-class classification using a multinomial generalized linear model (GLM).
#' It allows the user to define custom folds for cross-validation (e.g., leave-one-out if \code{folds = 1:nrow(datasc)}).
#' Optional SMOTE resampling can be applied to each training fold to correct for class imbalance.
#' Performance is assessed via confusion matrices and multiclass ROC curves.
#'
#' Requires that `class` is a factor and that `datasc` contains a `sample` column (which is removed internally).
#'
#' The ROC objects returned are computed using `pROC::multiclass.roc`, which provides an aggregate AUC estimate for multi-class settings.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  datasc <- my_data_frame
#'  levels <- c("A", "B", "C")
#'  vars <- c("gene1", "gene2", "gene3")
#'  result <- make_glm_l1o_multiclass(datasc, levs = levels, varnames = vars)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}},
#'  \code{\link[dplyr]{mutate}},
#'  \code{\link[assertthat]{assert_that}},
#'  \code{\link[pROC]{multiclass.roc}},
#'  \code{\link[nnet]{multinom}},
#'  \code{\link[UBL]{SmoteClassif}},
#'  \code{\link[caret]{confusionMatrix}}
#' @rdname make_glm_l1o_multiclass
#' @export
#' @importFrom dplyr select mutate
#' @importFrom assertthat assert_that
#' @importFrom pROC multiclass.roc
#' @importFrom caret confusionMatrix
#' @importFrom UBL SmoteClassif
#' @importFrom nnet multinom
make_glm_l1o_multiclass <- function(datasc, levs, varnames, folds=c(),
                                    do_smote=FALSE,
                                    smote_params=smote_params_default){
  predict_glm1 <- c()
  df <- datasc %>% dplyr::select(-sample) %>% dplyr::select(class, all_of(varnames))
  formula <- paste0("class ~ ", paste(varnames, sep="+", collapse="+")) %>% as.formula()

  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }

  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]

    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- UBL::SmoteClassif(form, smote_df,
                                C.perc = smote_params$dup_size,
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% dplyr::mutate(class=factor(class, levels=levs))
    }else{
      smoteData = NULL
    }

    mod_glm <- nnet::multinom(formula, data=train_df)
    newpred <- predict(mod_glm, test_df, type="prob")
    predict_glm1 <- rbind(predict_glm1, newpred)

  }

  predict1<- colnames(predict_glm1)[apply(predict_glm1, MAR=1, which.max)] %>% factor
  confmat1 <- confusionMatrix(predict1, factor(datasc$class))

  assertthat::assert_that(all(colnames(predict_glm1) == levels(datasc$class)))
  roc1 <- multiclass.roc(response=datasc$class, predictor=predict_glm1)
  roc_auc <- as.numeric(roc1$auc)

  mod_all <- nnet::multinom(formula, data=datasc, family = binomial)
  predict2 <- predict(mod_all, df, type="probs")
  classes2 <- colnames(predict2)[apply(predict2, MAR=1, \(x)which(x==max(x)))] %>% factor
  confmat2 <- confusionMatrix(classes2, factor(datasc$class))

  roc_obj_fullmod <- apply(predict2, MAR=2, \(x) multiclass.roc(datasc$class, x))
  roc_auc_fullmod <- sapply(roc_obj_fullmod, \(x)x$auc) %>% mean
  return(list(confmat=confmat1,
              confmat_no_l1o=confmat2,
              mod=mod_all,
              preds=predict1,
              preds_no_l1o=classes2,
              roc_obj_no_l1o=roc_obj_fullmod,
              roc_auc_no_l1o=roc_auc_fullmod,
              roc_obj=roc1,
              roc_auc=roc_auc
  ))
}
