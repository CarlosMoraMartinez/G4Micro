#' @title K-Nearest Neighbors Classification with Cross-Validation and Optional SMOTE
#' @description
#' Performs k-NN classification on a dataset with multiple values of k.
#' Supports cross-validation with any folds and optional SMOTE oversampling.
#' The input data must include a `sample` column.
#'
#' @param datasc A data frame containing the features, the class factor column named `class`, and a `sample` column.
#' @param levs A character vector of factor levels for the classification target (`class`).
#' @param varnames Character vector of variable names to be used as features.
#' @param different_ks Integer vector of k values to test for k-NN. Default: `c(1,3,5,7,9,11,13)`.
#' @param folds Integer vector of row indices to use as test folds in cross-validation. Default: all rows (leave-one-out).
#' @param do_smote Logical indicating whether to apply SMOTE oversampling on training data. Default: `FALSE`.
#' @param smote_params List with SMOTE parameters: `K` (neighbors) and `dup_size` (duplication strategy). Default: `list(K=5, dup_size="balance")`.
#'
#' @return A named list of results for each tested k, each containing:
#' \item{preds}{Predicted classes from cross-validation.}
#' \item{confmat}{Confusion matrix of cross-validated predictions.}
#' \item{pred_probs}{Predicted class probabilities.}
#' \item{preds_no_l1o}{Predictions from training on full data (no cross-validation).}
#' \item{confmat_no_l1o}{Confusion matrix of full-data predictions.}
#' \item{roc_obj}{ROC curve object.}
#' \item{roc_auc}{Area under the ROC curve.}
#'
#' @details
#' Cross-validation is done over the specified `folds` indices, which can represent any partitioning scheme.
#' SMOTE can be applied during training to handle class imbalance.
#' The function requires the packages `class`, `caret`, `pROC`, and optionally `smotefamily`.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(dplyr)
#'
#'   data(iris)
#'   # Add a sample column
#'   iris$sample <- 1:nrow(iris)
#'   # Convert Species to factor with levels
#'   iris$class <- factor(iris$Species)
#'   levels_iris <- levels(iris$class)
#'   # Use only first two classes for binary example:
#'   iris_bin <- iris %>% filter(class != "virginica") %>% droplevels()
#'   varnames_iris <- colnames(iris)[1:4]
#'   folds_iris <- sample(1:nrow(iris_bin), size = nrow(iris_bin), replace = FALSE)
#'
#'   res <- makeKnn_l1o(iris_bin, levs=levels(iris_bin$class), varnames=varnames_iris,
#'                      different_ks = c(3,5,7),
#'                      folds = folds_iris,
#'                      do_smote = FALSE)
#'   print(res$`K=3`$confmat)
#' }
#' }
#' @seealso
#' \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}, \code{\link[class]{knn}}, \code{\link[caret]{confusionMatrix}}, \code{\link[pROC]{roc}}
#' @rdname makeKnn_l1o
#' @export
#' @importFrom dplyr mutate select all_of
#' @importFrom caret confusionMatrix
#' @importFrom UBL SmoteClassif
#' @importFrom pROC roc multiclass.roc
#' @importFrom class knn
makeKnn_l1o <- function(datasc, levs, varnames,
                        different_ks=c(1, 3, 5, 7, 9, 11, 13),
                        folds=c(),
                        do_smote=FALSE,
                        smote_params=list(K=5, dup_size="balance")){
  results <- list()
  datasc <- datasc %>% dplyr::mutate(class=factor(class))
  train_df_all <- datasc %>%
    dplyr::select(-class, -sample) %>%
    dplyr::select(all_of(varnames))

  if(length(folds)==0){
    folds <- 1:nrow(datasc)
    reorder_samples <- FALSE
  }else{
    reorder_samples <- TRUE
  }
  for(k in  different_ks){
    kname = paste("K=", as.character(k), sep="", collapse="")
    results[[kname]] <- list()
    preds <- list()
    pred_probs <- list()
    for(i in folds){
      train_df <- train_df_all[-i, ]
      test_df <- train_df_all[i, ]

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

      kires <- class::knn(train_df, test_df, train_labels, k = k, prob = T)
      if(length(i) > 1){
        for(ii in 1:length(i)) preds[[ i[[ii]] ]] <- kires[ii]
        for(ii in 1:length(i)) pred_probs[[ i[[ii]] ]] <- class::knn(train_df, test_df[ii, ], train_labels, k = k, prob = T)
      }else{
        preds  <- c(preds, kires)
        pred_probs[[i]] <- kires
      }

    }# for i in folds

    preds <- unlist(preds)

    if(length(levs == 2)){
      prob_vec <- map_vec(pred_probs, \(x){ifelse(x==levs[1], 1-attr(x, "prob"), attr(x, "prob"))})
      roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=prob_vec)
      results[[kname]][["roc_obj"]] <- roc1
      results[[kname]][["roc_auc"]] <- as.numeric(roc1$auc)
    }else{
      prob_vec <- map_vec(pred_probs, \(x){ifelse(x==levs[1], 1-attr(x, "prob"), attr(x, "prob"))})
      roc1 <- multiclass.roc(response=datasc$class, predictor=prob_vec)
      results[[kname]][["roc_obj"]] <- roc1
      results[[kname]][["roc_auc"]] <- as.numeric(roc1$auc)
    }

    results[[kname]][["preds"]] <- levs[preds] %>% factor(levels=levs)
    results[[kname]][["confmat"]] <- confusionMatrix(results[[kname]][["preds"]],
                                                     factor(datasc$class),
                                                     positive = levs[2])
    results[[kname]][["pred_probs"]] <- prob_vec



    preds2 <- class::knn(train_df_all, train_df_all, datasc$class, k = k, prob = T)
    results[[kname]][["preds_no_l1o"]] <- levs[preds2] %>% factor(levels=levs)
    results[[kname]][["confmat_no_l1o"]] <- confusionMatrix(results[[kname]][["preds_no_l1o"]],
                                                            factor(datasc$class),
                                                            positive = levs[2])
  } # for k

  return(results)
}
