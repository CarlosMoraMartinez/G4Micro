#' @title Create Ensemble Model Predictions from Multiple Classifiers
#' @description This function performs ensemble classification using predictions
#' from multiple models, with optional performance-based weighting. It returns a confusion matrix,
#' predicted labels, and AUC for binary or multiclass classification.
#' @param datasc A data frame containing the true class labels in a column named `class`.
#' @param levs A character vector of class levels (e.g., \code{c("control", "case")}).
#' @param modlist A named list of trained models with prediction outputs in \code{modlist[[model]]$preds}.
#' @param model_res A data frame containing evaluation metrics for each model, including a column with model names (`model`) and a column with the performance metric specified in \code{param}.
#' @param param A character string indicating the name of the model performance column to use for model selection and weighting. Default is \code{"Kappa_l1out"}.
#' @param min_val Minimum acceptable value for the selected performance metric. Only models with at least this value will be used in the ensemble. Default: \code{0}.
#' @param prop Logical; if \code{TRUE}, models are weighted proportionally to their performance on \code{param}. If \code{FALSE}, all selected models have equal weight. Default: \code{TRUE}.
#' @param only_1_knn Logical; if \code{TRUE}, only the best KNN model is kept (others are removed). Default: \code{FALSE}.
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{confmat}: Confusion matrix from the ensemble predictions.
#'   \item \code{confmat_no_l1o}: Always \code{NULL} (placeholder).
#'   \item \code{mod}: Always \code{NULL} (placeholder).
#'   \item \code{preds}: Factor vector of predicted class labels.
#'   \item \code{pred_df}: Data frame of individual model predictions used for the ensemble.
#'   \item \code{preds_no_l1o}: Always \code{NULL} (placeholder).
#'   \item \code{roc_obj_no_l1o}: Always \code{NULL} (placeholder).
#'   \item \code{roc_auc_no_l1o}: Always \code{NULL} (placeholder).
#'   \item \code{roc_obj}: ROC object from \code{pROC::roc()} or \code{pROC::multiclass.roc()}.
#'   \item \code{roc_auc}: Numeric AUC value from the ensemble prediction.
#' }
#' @details
#' The function aggregates predictions from multiple models based on either equal weighting or performance-weighted voting.
#' If \code{only_1_knn = TRUE}, only the top KNN model (based on \code{param}) is used; others are excluded.
#'
#' In binary classification, the ensemble probability is used to compute AUC with \code{pROC::roc()}. For multiclass problems, \code{pROC::multiclass.roc()} is used.
#'
#' All models must have predictions available in \code{modlist[[model]]$preds}. The true class labels must be in \code{datasc$class}.
#'
#' Typical values for \code{param} include:
#' \itemize{
#'   \item \code{"Kappa_l1out"}: Cohen's kappa on CV predictions.
#'   \item \code{"Accuracy_l1out"}: Accuracy on CV predictions.
#'   \item \code{"AUC_l1out"}: AUC on CV predictions.
#' }
#'
#' Models are selected if their \code{param} is greater than or equal to \code{min_val}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{filter}}
#'  \code{\link[purrr]{map}}
#' @rdname make_ensemble_votes
#' @export
#' @importFrom dplyr arrange filter
#' @importFrom purrr map
#' @importFrom pROC roc multiclass.roc
#' @importFrom caret confusionMatrix
make_ensemble_votes <- function(datasc, levs, modlist, model_res, param="Kappa_l1out",
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
  preddf <- map(m2use, \(x) tibble( !!x := modlist[[x]]$preds))  %>% bind_cols()
  #names(preddf) <- m2use
  if(prop){
    ponderfac <- model_res[match(m2use, model_res$model), param]
    ponderfac <- (ponderfac - min(ponderfac))/(max(ponderfac) - min(ponderfac)) + 0.1
  }else{
    ponderfac <- rep(1, length(m2use))
  }
  preds <- c()
  votes <- list()
  for(i in 1:nrow(preddf)){
    classcore <- map_vec(levs, \(ll) sum(ponderfac[preddf[i, ] == ll]))
    names(classcore) <- levs
    l1 <- length(preds)
    preds <- c(preds, levs[which.max(classcore)] )
    l2 <- length(preds)
    cat(i, ": L1=", l1, ", L2=", l2, ifelse(l1==l2, " --WARNING--", ""),  "\n")
    votes[[i]] <- classcore
  }
  preds <- factor(preds, levels=levs)
  confmat1 <- confusionMatrix(preds, datasc$class, positive = levs[2])

  if(length(levs) == 2){
    probs <- map_vec(votes, \(x)x[levs[2]]/sum(x) )
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs)
  }else{
    probs <- purrr::map(votes, .f = \(x) x/sum(x)) %>%
      bind_rows %>% as.matrix #%>% pull(!!sym(levs[2]))
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs)
  }
  return(list(confmat=confmat1,
              confmat_no_l1o=NULL,
              mod=NULL,
              preds=preds,
              pred_df=preddf,
              preds_no_l1o=NULL,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}
