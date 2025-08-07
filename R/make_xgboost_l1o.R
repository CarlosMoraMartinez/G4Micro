#' @title Train and Evaluate XGBoost Model with Cross-Validation and Optional SMOTE Oversampling
#' @description
#' This function trains an XGBoost classification model using generalised cross-validation.
#' It supports optional class balancing via weights and optional SMOTE oversampling.
#' The function returns performance metrics including confusion matrices and ROC curves.
#'
#' @param datasc A data frame containing the dataset with predictors, a `class` factor column, and optionally a `sample` column.
#' @param levs A character vector of factor levels for the target classification variable `class`.
#' @param varnames A character vector specifying the names of predictor variables to be used.
#' @param folds Integer vector specifying which samples to leave out for L1O; default uses all samples (leave-one-out).
#' @param xgboost_params A named list of parameters to control XGBoost model training (e.g., max_depth, learning_rate, nrounds, etc.). Default: \code{list(
#'   learning_rate = 0.3,
#'   max_depth = 2,
#'   nrounds = 30,
#'   min_child_weight = 1,
#'   subsample = 1,
#'   colsample_bytree = 0.6,
#'   reg_lambda = 1,
#'   gamma = 0,
#'   reg_alpha = 0,
#'   nthread = 1,
#'   objective = "binary:logistic",
#'   balance_weights = TRUE
#' )}
#' @param do_smote Logical indicating whether to perform SMOTE oversampling on training data folds; default is FALSE.
#' @param smote_params A named list of parameters for SMOTE (e.g., K = 5, dup_size = "balance").
#'
#' @return A list containing:
#' \item{confmat}{Confusion matrix from CV predictions}
#' \item{confmat_no_l1o}{Confusion matrix from model trained on full dataset}
#' \item{mod}{XGBoost model trained on full dataset}
#' \item{preds}{Predicted classes from CV}
#' \item{pred_probs}{Predicted probabilities from CV}
#' \item{preds_no_l1o}{Predicted classes from model trained on full dataset}
#' \item{roc_obj_no_l1o}{ROC object from full dataset model predictions (currently NULL)}
#' \item{roc_auc_no_l1o}{ROC AUC from full dataset model (currently NULL)}
#' \item{roc_obj}{ROC object from CV predictions}
#' \item{roc_auc}{ROC AUC from CV predictions}
#' \item{xgboost_params}{Parameters used for XGBoost training}
#' \item{smoteData}{Last SMOTE dataset generated (or NULL if not used)}
#'
#' @details
#' - This function uses leave-one-out cross-validation by default (each sample is left out once).
#' - SMOTE can be enabled for oversampling the minority class in training folds.
#' - Class balancing weights are automatically calculated if enabled.
#' - Currently, ROC objects for the model trained on the full dataset are not computed.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   data(iris)
#'   iris <- iris[iris$Species != "setosa", ]
#'   iris$class <- factor(iris$Species)
#'   iris$sample <- 1:nrow(iris)
#'   params <- list(max_depth = 3, learning_rate = 0.1, nrounds = 50,
#'                  min_child_weight = 1, subsample = 0.8, colsample_bytree = 0.8,
#'                  gamma = 0, reg_lambda = 1, reg_alpha = 0, nthread = 1,
#'                  objective = "binary:logistic", balance_weights = TRUE)
#'   result <- make_xgboost_l1o(iris, levs = levels(iris$class),
#'                              varnames = names(iris)[1:4],
#'                              xgboost_params = params,
#'                              do_smote = FALSE)
#'   print(result$confmat)
#' }
#' }
#'
#' @seealso
#'  \code{\link[dplyr]{select}},
#'  \code{\link[xgboost]{xgboost}},
#'  \code{\link[pROC]{roc}},
#'  \code{\link[caret]{confusionMatrix}},
#'  \code{\link[UBL]{SmoteClassif}},
#'
#' @rdname make_xgboost_l1o
#' @export
#' @importFrom dplyr select
#' @importFrom xgboost xgboost
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc multiclass.roc
#' @importFrom UBL SmoteClassif
#'
make_xgboost_l1o <- function(datasc, levs, varnames,
                             folds=folds(),
                             xgboost_params = xgboost_params_default,
                             do_smote=FALSE,
                             smote_params=smote_params_default
){
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  predict_probs = numeric(0)

  if(xgboost_params$balance_weights){

    #class_weights <- class_weights/min(class_weights)
    #weights <- class_weights[datasc$class]

    class_weights <- 1/table(datasc$class)
    if(length(levs)==2){
      posweight <- class_weights[levs[1]]/class_weights[levs != levs[1]]
    } else{
      posweight <- as.vector(class_weights)
      names(posweight) <- names(class_weights)
    }

  }else{
    class_weights <- NULL
    if(length(levs)==2){
      posweight <- 1
    }else{
      posweight <- rep(1, length(levs))
      names(posweight) <- levs
    }
  }

  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]

    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    #if(xgboost_params$balance_weights){
    #  train_weights <- weights[-i]
    #}else{
    #  train_weights <- NULL
    #}

    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df,
                                C.perc = smote_params$dup_size,
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
      #if(xgboost_params$balance_weights){
      #  train_weights <- class_weights[train_labels]
      #}
    }else{
      smoteData = NULL
    }
    mod_tree1 <- xgboost(x = train_df, y = train_labels,
                         #weights=train_weights,
                         scale_pos_weight = posweight,
                         max_depth = xgboost_params$max_depth,
                         learning_rate = xgboost_params$learning_rate,
                         nrounds = xgboost_params$nrounds,
                         min_child_weight = xgboost_params$min_child_weight,
                         subsample = xgboost_params$subsample,
                         colsample_bytree = xgboost_params$colsample_bytree,
                         gamma = xgboost_params$gamma,
                         reg_lambda=xgboost_params$reg_lambda,
                         reg_alpha= xgboost_params$reg_alpha,
                         nthread = xgboost_params$nthread,
                         objective = xgboost_params$objective)
    if(length(levs) == 2){
      predict_probs <- c(predict_probs, predict(mod_tree1, test_df))
    }else{
      predict_probs <- rbind(predict_probs, predict(mod_tree1, test_df))
    }


  }
  if(length(levs) == 2){
    predict_tree1 <- factor(levs[as.integer(round(predict_probs))+1], levels=levs)
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=predict_probs)
  }else{
    predict_tree1 <- factor(colnames(predict_probs)[apply(predict_probs, MAR=1, which.max)], levels=levs)
    roc1 <- multiclass.roc(response=datasc$class, predictor=predict_probs)
  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])

  mod_tree1 <-  xgboost(x = df, y = datasc$class,
                        #weights=weights,
                        scale_pos_weight = posweight,
                        max_depth = xgboost_params$max_depth,
                        learning_rate = xgboost_params$learning_rate,
                        nrounds = xgboost_params$nrounds,
                        min_child_weight = xgboost_params$min_child_weight,
                        subsample = xgboost_params$subsample,
                        colsample_bytree = xgboost_params$colsample_bytree,
                        gamma = xgboost_params$gamma,
                        reg_lambda=xgboost_params$reg_lambda,
                        reg_alpha= xgboost_params$reg_alpha,
                        nthread = xgboost_params$nthread,
                        objective = xgboost_params$objective)

  predict_tree2 <- predict(mod_tree1, df)
  if(length(levs) == 2){
    predict_tree2 <- factor(levs[as.integer(round(predict_tree2))+1], levels=levs)
  }else{
    predict_tree2 <- factor(levs[apply(predict_tree2, MAR=1, which.max)], levels=levs)
  }
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])

  return(list(confmat=confmat_tree1,
              confmat_no_l1o=confmat_tree2,
              mod=mod_tree1,
              preds=predict_tree1,
              pred_probs = predict_probs,
              preds_no_l1o=predict_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc),
              xgboost_params = xgboost_params,
              smoteData=smoteData))
}
