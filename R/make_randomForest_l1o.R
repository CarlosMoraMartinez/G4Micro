#' @title Train and Evaluate Random Forest with Leave-One-Out Cross-Validation
#' @description
#' This function trains a Random Forest classifier on the provided dataset using
#' leave-one-out or k-fold cross-validation. It optionally applies SMOTE for class balancing
#' during training folds. The function returns performance metrics including confusion matrices,
#' ROC curves, and the trained model on the full dataset.
#'
#' @param datasc A data frame containing the dataset with features and a factor column named `class` indicating class labels, and a column named `sample`.
#' @param levs A character vector of factor levels for the class variable. The second level is considered the "positive" class.
#' @param varnames A character vector specifying which columns in `datasc` to use as predictor variables.
#' @param folds A vector of indices indicating which samples to leave out in each iteration. If empty, leave-one-out is used. Default: `c()`.
#' @param randomforest_params A list of parameters for the randomForest model, including `ntree`, `mtry`, `nodesize`, and `balance_weights` (logical).
#'                            Default is `randomforest_params_default`.
#' @param do_smote Logical, whether to apply SMOTE oversampling on the training folds to balance classes. Default is FALSE.
#' @param smote_params A list of parameters for SMOTE, including `K` (number of neighbors) and `dup_size` (duplication size or "balance"). Default is `smote_params_default`.
#'
#' @return A list containing:
#' \describe{
#'   \item{confmat}{Confusion matrix for the cross-validation predictions.}
#'   \item{confmat_no_l1o}{Confusion matrix for predictions from the model trained on the entire dataset.}
#'   \item{mod}{Random Forest model trained on the entire dataset.}
#'   \item{preds}{Vector of predicted classes from CV.}
#'   \item{pred_probs}{Numeric vector or data frame of predicted probabilities from CV.}
#'   \item{preds_no_l1o}{Predicted classes from the full model.}
#'   \item{roc_obj}{ROC curve object from CV predictions.}
#'   \item{roc_auc}{AUC value from CV predictions.}
#'   \item{smoteData}{The last SMOTE-augmented training set generated during CV (if `do_smote = TRUE`).}
#'   \item{params}{List of parameters used for the Random Forest model.}
#' }
#'
#' @details
#' The function supports binary and multiclass classification. If `balance_weights` is TRUE in `randomforest_params` and `do_smote` is FALSE,
#' class weights are set inversely proportional to class frequencies.
#' If `do_smote` is TRUE, SMOTE oversampling is applied on each training fold with parameters specified in `smote_params`.
#' The function requires the `randomForest` package and optionally the `UBL` package for SMOTE.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example usage:
#'   data(iris)
#'   iris <- iris %>%
#'   mutate(class = factor(ifelse(Species == "setosa", "pos", "neg")), sample = 1:n())
#'
#'   rf_params <- list(ntree=500, mtry=3, nodesize=5, balance_weights=TRUE)
#'   smote_params <- list(K=5, dup_size="balance")
#'
#'   result <- make_randomForest_l1o(datasc=iris, levs=c("neg", "pos"),
#'                      varnames=colnames(iris)[1:4], folds=1:nrow(iris),
#'                      randomforest_params=rf_params, do_smote=TRUE,
#'                      smote_params=smote_params)
#'   print(result$confmat)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{select}},
#' \code{\link[randomForest]{randomForest}},
#' \code{\link[UBL]{SMOTE}}
#'
#' @rdname make_randomForest_l1o
#' @export
#' @importFrom dplyr select
#' @importFrom UBL SmoteClassif
#' @importFrom randomForest randomForest
#' @importFrom pROC roc multiclass.roc
#' @importFrom caret confusionMatrix
make_randomForest_l1o <- function(datasc, levs, varnames,
                                  folds=folds(),
                                  randomforest_params = randomforest_params_default,
                                  do_smote=FALSE,
                                  smote_params=smote_params_default
){
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))

  if(randomforest_params$balance_weights & !do_smote){
    classweights <- table(datasc$class)
    sweights <- max(classweights)/classweights[datasc$class]
  }else{
    sweights <- rep(1, nrow(datasc))
  }


  predict_tree1 <- factor()
  predict_tree1_probs <- list()


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
    train_weighs <- sweights[-i]

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
      train_weighs <- rep(1, nrow(train_df))
    }else{
      smoteData = NULL
    }
    mod_tree1 <- randomForest(x=train_df, y=train_labels, levels=levs,
                              weights = train_weighs,
                              ntree = randomforest_params$ntree,
                              mtry = randomforest_params$mtry,
                              nodesize = randomforest_params$nodesize)
    predict_tree1 <- c(predict_tree1, predict(mod_tree1, test_df))
    if(length(i) > 1){
      for(ii in i) predict_tree1_probs[[ii]] <- predict(mod_tree1, test_df[as.character(ii), ], type = "prob")
    }else{
      predict_tree1_probs[[i]] <- predict(mod_tree1, test_df, type = "prob")
    }

  }
  if(reorder_samples){
    predict_tree1 <- predict_tree1[as.character(1:nrow(datasc))]
  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict_tree1_probs, \(xx) xx %>% as.data.frame) %>%
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- predict_tree1_probs %>% map(as.data.frame) %>% bind_rows
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
    roc_auc <- as.numeric(roc1$auc)
  }

  mod_tree1 <- randomForest(x=df, y=datasc$class, levels=levs,
                            weights = sweights,
                            ntree = randomforest_params$ntree,
                            mtry = randomforest_params$mtry,
                            nodesize = randomforest_params$nodesize)
  predict_tree2 <- predict(mod_tree1, df)
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  return(list(confmat=confmat_tree1,
              confmat_no_l1o=confmat_tree2,
              mod=mod_tree1,
              preds=predict_tree1,
              pred_probs =probs_vector,
              preds_no_l1o=predict_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc),
              smoteData=smoteData,
              params =  randomforest_params))
}
