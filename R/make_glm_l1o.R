#' @title Generalized GLM Cross-Validation with Optional SMOTE
#' @description
#' Performs logistic regression with cross-validation based on user-defined folds,
#' optionally using SMOTE to address class imbalance in the training data.
#' Returns predictions, performance metrics, and ROC curves for both the cross-validated and the full-model fits.
#'
#' @param datasc A data frame containing the `sample`, `class`, and feature variables.
#' @param levs A character vector with two levels defining the factor levels for the classification task.
#' @param varnames A character vector with the names of predictor variables to use in the model.
#' @param folds A numeric vector specifying row indices for the test set in each iteration of cross-validation. Defaults to leave-one-out: `1:nrow(datasc)`.
#' @param do_smote Logical, whether to apply SMOTE to the training data. Default: `FALSE`.
#' @param smote_params A list with SMOTE parameters: `K` (number of neighbors) and `dup_size` (oversampling rate). Default: `list(K = 5, dup_size = 2)`.
#'
#' @return A list containing confusion matrices, predictions, prediction probabilities, ROC objects, and AUC values for both the cross-validated and full models.
#'
#' @details
#' This function builds a logistic regression model using the specified variables and cross-validation strategy.
#' If `do_smote = TRUE`, SMOTE is applied to the training set at each iteration.
#' The function returns metrics such as accuracy, sensitivity, specificity, and AUC for both the cross-validated predictions and the model trained on the full dataset.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   res <- make_glm_l1o(datasc, levs = c("Control", "Case"),
#'                      varnames = c("feature1", "feature2"),
#'                      folds = 1:10,
#'                      do_smote = TRUE)
#'   print(res$confmat)
#' }
#' }
#'
#' @seealso
#'  \code{\link[dplyr]{select}},
#'  \code{\link[dplyr]{mutate}},
#'  \code{\link[UBL]{SmoteClassif}},
#'  \code{\link[stats]{glm}},
#'  \code{\link[stats]{predict}},
#'  \code{\link[pROC]{roc}},
#'  \code{\link[caret]{confusionMatrix}}
#'
#' @rdname make_glm_l1o
#' @export
#' @importFrom dplyr select mutate
#' @importFrom UBL SmoteClassif
#' @importFrom stats glm predict as.formula
#' @importFrom pROC roc
#' @importFrom caret confusionMatrix
make_glm_l1o <- function(datasc, levs, varnames, folds= c(),
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
      smoteData <- SmoteClassif(form, smote_df,
                                C.perc = smote_params$dup_size,
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% dplyr::mutate(class=factor(class, levels=levs))
    }else{
      smoteData = NULL
    }

    mod_glm <- glm(formula, data=train_df, family = binomial)
    predict_glm1 <- c(predict_glm1, predict(mod_glm, test_df, type = "response"))

  }
  predict1 <- ifelse(predict_glm1 > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
  confmat1 <- confusionMatrix(predict1, datasc$class, positive = levs[2])
  roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=predict_glm1)

  mod_all <- glm(formula, data=datasc, family = binomial)
  predict2_probs <- predict(mod_all, df, type = "response")
  predict2 <- ifelse(predict2_probs > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
  confmat2 <- confusionMatrix(predict2, datasc$class, positive = levs[2])
  roc2 <- roc(response=as.numeric(datasc$class)-1, predictor=predict2_probs)
  return(list(confmat=confmat1,
              confmat_no_l1o=confmat2,
              mod=mod_all,
              preds=predict1,
              pred_probs = predict_glm1,
              preds_no_l1o=predict2,
              pred_probs_no_l1o=predict2_probs,
              roc_obj_no_l1o=roc2,
              roc_auc_no_l1o=as.numeric(roc2$auc),
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}
