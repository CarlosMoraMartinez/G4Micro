#' @title Train CatBoost model with cross-validation and optional SMOTE augmentation
#' @description
#' Performs cross-validation training of a CatBoost classifier with support for categorical imbalance handling,
#' optional SMOTE oversampling on training folds, and returns predictions, probabilities, confusion matrices,
#' and trained models. Supports binary and multiclass classification.
#'
#' @param datasc Data frame containing the dataset. Must include a factor column named `class` and a column `sample`.
#' @param levs Character vector specifying the levels of the class factor, defining class order.
#' @param varnames Character vector of variable names (features) to use as predictors.
#' @param folds Integer vector specifying indices for cross-validation folds. Default is all rows (no CV).
#' @param catboost_params A named list of CatBoost training parameters controlling model training behavior.
#'   Defaults to:
#'   \itemize{
#'     \item \code{iterations} (int): Number of boosting iterations. Default: 100
#'     \item \code{learning_rate} (numeric): Step size shrinkage used to prevent overfitting. Default: 0.05
#'     \item \code{depth} (int): Depth of the tree. Default: 2
#'     \item \code{loss_function} (string): Loss function to optimize. Default: "Logloss"
#'     \item \code{eval_metric} (string): Metric for evaluation. Default: "AUC"
#'     \item \code{random_seed} (int): Seed for reproducibility. Default: 123
#'     \item \code{use_best_model} (bool): Use best model found during training. Default: TRUE
#'     \item \code{od_type} (string): Overfitting detector type. Default: "Iter"
#'     \item \code{od_wait} (int): Overfitting detector wait iterations. Default: 20
#'     \item \code{verbose} (bool): Print training progress. Default: FALSE
#'     \item \code{thread_count} (int): Number of threads to use. Default: 1
#'     \item \code{balance_weights} (bool): Whether to balance class weights. Default: TRUE
#'     \item \code{bootstrap_type} (string): Bootstrap sampling type ("Bayesian", "Bernoulli", etc.). Default: "Bayesian"
#'     \item \code{l2_leaf_reg} (numeric): L2 regularization coefficient. Default: 3
#'     \item \code{subsample} (numeric): Subsample ratio for Bernoulli bootstrap (not used if "Bayesian" is selected). Default: 0.6
#'     \item \code{grow_policy} (string): Tree grow policy ("Depthwise" or others). Default: "Depthwise"
#'     \item \code{auto_class_weights} (string or bool): Automatic class weights adjustment. Default: "Balanced"
#'   }
#' @param do_smote Logical indicating whether to apply SMOTE oversampling on training folds. Default is FALSE.
#' @param smote_params List with parameters for SMOTE sampling, including `K` (nearest neighbors) and `dup_size` (duplication size or balance strategy). Default is list(K = 5, dup_size = "balance").
#'
#' @return A list containing:
#' \item{confmat}{Confusion matrix of predictions from cross-validation folds.}
#' \item{confmat_no_l1o}{Confusion matrix of predictions from model trained on full data (no CV).}
#' \item{mod}{CatBoost model trained on the full dataset.}
#' \item{preds}{Factor vector of predicted classes from cross-validation folds.}
#' \item{pred_probs}{Numeric vector or matrix of predicted class probabilities from cross-validation folds.}
#' \item{preds_no_l1o}{Factor vector of predicted classes from the full model.}
#' \item{roc_obj_no_l1o}{Currently NULL; placeholder for ROC object from full model predictions.}
#' \item{roc_auc_no_l1o}{Currently NULL; placeholder for AUC of full model predictions.}
#' \item{roc_obj}{ROC object from cross-validation predictions (binary or multiclass).}
#' \item{roc_auc}{Numeric AUC value from cross-validation predictions.}
#' \item{xgboost_params}{Parameters used for CatBoost training.}
#' \item{smoteData}{The last SMOTE-augmented training data frame (if `do_smote = TRUE`), otherwise NULL.}
#'
#' @details
#' The function trains a CatBoost classification model using the specified cross-validation folds.
#' It supports binary and multiclass classification by encoding labels appropriately.
#' If `do_smote` is TRUE, SMOTE oversampling is applied to training folds to balance classes before model training.
#' Class weights can be used if specified in `catboost_params` and SMOTE is not applied.
#' After cross-validation, a final CatBoost model is trained on the entire dataset without folds.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   data(iris)
#'   df <- %>% filter(Species != "setosa")
#'   df$class <- df$Species
#'   df$sample <- 1:nrow(df)
#'   levs <- as.character(unique(df$class))
#'   varnames <- colnames(df)[1:4]
#'   folds <- caret::createFolds(df$class, k = 10, list = TRUE, returnTrain = FALSE) # 10-fold CV
#'   params <- list(depth=6,
#'                  learning_rate=0.05,
#'                  iterations=100,
#'                  loss_function="Logloss",
#'                  eval_metric="AUC",
#'                  bootstrap_type="Bernoulli",
#'                  random_seed = 123,
#'                  use_best_model = TRUE,
#'                  od_type = "Iter",
#'                  od_wait = 20,
#'                  l2_leaf_reg=3,
#'                  subsample=0.8,
#'                  grow_policy="Depthwise",
#'                  auto_class_weights="Balanced",
#'                  verbose=FALSE,
#'                  thread_count=4,
#'                  balance_weights=TRUE)
#'   results <- make_catboost_l1o(df, levs, varnames, folds=folds, catboost_params=params, do_smote=TRUE)
#'   print(results$confmat)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}
#'  \code{\link[catboost]{catboost.train}}
#'  \code{\link[pROC]{roc}},
#'  \code{\link[caret]{confusionMatrix}},
#'  \code{\link[UBL]{SmoteClassif}},
#' @rdname make_catboost_l1o
#' @export
#' @importFrom dplyr select
#' @importFrom catboost catboost.load_pool catboost.train catboost.predict
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc multiclass.roc
#' @importFrom UBL SmoteClassif
make_catboost_l1o <- function(datasc, levs, varnames,
                              folds=folds(),
                              catboost_params = catboost_params,
                              do_smote=FALSE,
                              smote_params=list(K=5, dup_size="balance")
){
  datasc$class <- factor(datasc$class, levels=levs)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
    reorder_samples <- FALSE
  }else{
    reorder_samples <- TRUE
  }
  predict_probs = list()

  if(catboost_params$balance_weights & ! do_smote){
    class_weights <- table(datasc$class)
    class_weights_vec <- max(class_weights)/class_weights %>% as.vector
    names(class_weights_vec) <- levs
  } else {
    class_weights_vec <- rep(1, length(levs))
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
      #smoteData <- SMOTE(train_df, train_labels, K=smote_params$K, dup_size = smote_params$dup_size)
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
    if(length(levs) == 2){
      train_pool <- catboost::catboost.load_pool(data = train_df, label = as.integer(train_labels == levs[2]))
    }else{
      train_pool <- catboost::catboost.load_pool(data = train_df, label = as.integer(train_labels)-1)
    }
    test_pool <- catboost::catboost.load_pool(data = test_df)

    if(catboost_params$bootstrap_type == "Bernoulli"){
      model <- catboost::catboost.train(learn_pool = train_pool, params = list(
        depth = catboost_params$depth,
        learning_rate = catboost_params$learning_rate,
        iterations = catboost_params$iterations,
        loss_function = catboost_params$loss_function,
        eval_metric = catboost_params$eval_metric,
        bootstrap_type = catboost_params$bootstrap_type,
        l2_leaf_reg = catboost_params$l2_leaf_reg,
        subsample = catboost_params$subsample,
        grow_policy = catboost_params$grow_policy,
        auto_class_weights= catboost_params$auto_class_weights,
        thread_count = catboost_params$thread_count,
        #class_weights = class_weights_vec,
        logging_level = "Silent"
      ))
    }else{
      model <- catboost::catboost.train(learn_pool = train_pool, params = list(
        depth = catboost_params$depth,
        learning_rate = catboost_params$learning_rate,
        iterations = catboost_params$iterations,
        loss_function = catboost_params$loss_function,
        eval_metric = catboost_params$eval_metric,
        bootstrap_type = catboost_params$bootstrap_type,
        l2_leaf_reg = catboost_params$l2_leaf_reg,
        #subsample = catboost_params$subsample,
        grow_policy = catboost_params$grow_policy,
        auto_class_weights= catboost_params$auto_class_weights,
        thread_count = catboost_params$thread_count,
        #class_weights = class_weights_vec,
        logging_level = "Silent"
      ))
    }
    pred_prob <- catboost::catboost.predict(model, test_pool, prediction_type = "Probability")

    if(reorder_samples){
      if(length(levs) == 2){
        for(ii in 1:length(i)) predict_probs[[ i[ii] ]] <- pred_prob[ii]
      }else{
        rownames(pred_prob) <- i
        predict_probs <- rbind(predict_probs, pred_prob)
      }
    }else{
      if(length(levs) == 2){
        predict_probs <- c(predict_probs, pred_prob)
      }else{
        predict_probs <- rbind(predict_probs, pred_prob)
      }
    }

  }
  if(length(levs) == 2){
    predict_probs <- unlist(predict_probs)
  }else if(reorder_samples & length(levs) >2){
    predict_probs <- predict_probs[as.character(1:nrow(datasc)),] ## test late
  }


  if(length(levs) == 2){
    predict_tree1 <- factor(levs[as.integer(round(predict_probs))+1], levels=levs)
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=predict_probs)
  }else{
    predict_tree1 <- factor(levs[apply(predict_probs, MAR=1, which.max)], levels=levs)
    colnames(predict_probs) <- as.character(levels(datasc$class))
    roc1 <- multiclass.roc(response=datasc$class, predictor=predict_probs)
  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])

  train_pool <- catboost::catboost.load_pool(data = df, label = as.integer(datasc$class == levs[2]))
  test_pool <- catboost::catboost.load_pool(data = df)

  if(catboost_params$bootstrap_type == "Bernoulli"){
    model2 <- catboost::catboost.train(learn_pool = train_pool, params = list(
      depth = catboost_params$depth,
      learning_rate = catboost_params$learning_rate,
      iterations = catboost_params$iterations,
      loss_function = catboost_params$loss_function,
      eval_metric = catboost_params$eval_metric,
      bootstrap_type = catboost_params$bootstrap_type,
      l2_leaf_reg = catboost_params$l2_leaf_reg,
      subsample = catboost_params$subsample,
      grow_policy = catboost_params$grow_policy,
      auto_class_weights= catboost_params$auto_class_weights,
      thread_count = catboost_params$thread_count,
      #class_weights = class_weights_list,
      logging_level = "Silent"
    ))
  }else{
    model2 <- catboost::catboost.train(learn_pool = train_pool, params = list(
      depth = catboost_params$depth,
      learning_rate = catboost_params$learning_rate,
      iterations = catboost_params$iterations,
      loss_function = catboost_params$loss_function,
      eval_metric = catboost_params$eval_metric,
      bootstrap_type = catboost_params$bootstrap_type,
      l2_leaf_reg = catboost_params$l2_leaf_reg,
      #subsample = catboost_params$subsample,
      grow_policy = catboost_params$grow_policy,
      auto_class_weights= catboost_params$auto_class_weights,
      thread_count = catboost_params$thread_count,
      #class_weights = class_weights_list,
      logging_level = "Silent"
    ))
  }
  predict_tree2 <- catboost::catboost.predict(model2, test_pool, prediction_type = "Probability")
  if(length(levs) == 2){
    predict_tree2 <- factor(levs[as.integer(round(predict_tree2))+1], levels=levs)
  }else{
    predict_tree2 <- factor(levs[apply(predict_tree2, MAR=1, which.max)], levels=levs)
  }
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])

  return(list(confmat=confmat_tree1,
              confmat_no_l1o=confmat_tree2,
              mod=model2,
              preds=predict_tree1,
              pred_probs = predict_probs,
              preds_no_l1o=predict_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc),
              xgboost_params = catboost_params,
              smoteData=smoteData))
}
