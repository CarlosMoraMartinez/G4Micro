#' @title Combine Probabilistic Predictions from Multiple Models Using Cross-Validated Performance
#' @description This function combines probabilistic outputs from several classification models into
#' a single prediction using weighted averaging. Model performance is used as a weight, and predictions
#' are evaluated using cross-validation results.
#' @param datasc A data frame containing the samples used for prediction. Must include the actual class labels in a column named `class`.
#' @param levs A character vector of class levels (e.g., \code{c("control", "case")}).
#' @param modlist A named list of trained models, where each element is a list containing \code{$pred_probs}, a numeric vector
#' of predicted probabilities (typically from cross-validation).
#' @param model_res A data frame with model performance metrics from cross-validation. Must include a column named \code{model} and
#' a metric column used in \code{param}.
#' @param param A string with the name of the performance metric column in \code{model_res} to be used as
#' weight (e.g., 'BalancedAccuracy_l1out'). Default: \code{"BalancedAccuracy_l1out"}
#' @param min_val Minimum value of the chosen performance metric for a model to be included in the ensemble. Default: \code{0}
#' @param prop Logical; whether to weight predictions proportionally to the model's performance. Default: \code{TRUE}
#' @param only_1_knn Logical; if \code{TRUE}, only the best KNN model is used in the ensemble. Default: \code{FALSE}
#' @return A list with elements:
#' \itemize{
#'   \item \code{confmat}: Confusion matrix comparing predicted vs actual classes.
#'   \item \code{confmat_no_l1o}: Placeholder for compatibility (NULL).
#'   \item \code{mod}: Placeholder for compatibility (NULL).
#'   \item \code{preds}: Factor vector of predicted class labels.
#'   \item \code{pred_probs}: Numeric vector of averaged predicted probabilities.
#'   \item \code{pred_df}: Data frame of individual model probabilities.
#'   \item \code{preds_no_l1o}: Placeholder for compatibility (NULL).
#'   \item \code{roc_obj_no_l1o}: Placeholder for compatibility (NULL).
#'   \item \code{roc_auc_no_l1o}: Placeholder for compatibility (NULL).
#'   \item \code{roc_obj}: ROC object computed from ensemble probabilities.
#'   \item \code{roc_auc}: AUC value of the ensemble prediction.
#' }
#' @details
#' Models are selected based on their cross-validation performance (e.g., Balanced Accuracy),
#' and probabilistic predictions are combined using a weighted average. If \code{prop = TRUE},
#' weights are scaled to a 0–1 range and slightly regularized (+0.1).
#'
#' The ensemble is evaluated using the ROC AUC and confusion matrix, assuming binary classification.
#'
#' Note: This function assumes that \code{pred_probs} are numeric probabilities between 0 and 1, and class labels are binary and ordered.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{filter}}
#' @rdname make_ensemble_probs
#' @export
#' @importFrom dplyr arrange filter
#' @importFrom pROC roc multiclass.roc
#' @importFrom caret confusionMatrix
make_ensemble_probs <- function(datasc, levs, modlist, model_res, param="BalancedAccuracy_l1out", #"Kappa_l1out",
                                min_val=0, prop=TRUE,
                                only_1_knn=FALSE){ # 0.65
  remove_knn <- model_res %>%
    dplyr::arrange(desc(!!sym(param))) %>%
    dplyr::filter(grepl("KNN", model)) %>% pull(model)
  remove_knn <- remove_knn[2:length(remove_knn)]
  m2use <- model_res %>%
    dplyr::filter(!!sym(param) >= min_val) %>%
    dplyr::filter(! (model %in% remove_knn & only_1_knn)) %>%
    pull(model)
  preddf <- map(m2use, \(x) modlist[[x]]$pred_probs)  %>% bind_cols()
  names(preddf) <- m2use
  if(prop){
    ponderfac <- model_res[match(m2use, model_res$model), param]
    ponderfac <- (ponderfac - min(ponderfac))/(max(ponderfac) - min(ponderfac)) + 0.1

  }else{
    ponderfac <- rep(1, length(m2use))
  }
  preds <- c()
  avg_probs <- c()
  for(i in 1:nrow(preddf)){
    classcore <-sum(preddf[i, ]*ponderfac)/sum(ponderfac)
    avg_probs <- c(avg_probs, classcore)
    preds <- c(preds, levs[as.integer(round(classcore))+1] )
  }
  preds <- factor(preds, levels=levs)
  confmat1 <- confusionMatrix(preds, datasc$class, positive = levs[2])

  roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=avg_probs)

  return(list(confmat=confmat1,
              confmat_no_l1o=NULL,
              mod=NULL,
              preds=preds,
              pred_probs=avg_probs,
              pred_df=preddf,
              preds_no_l1o=NULL,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}
