#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param levs PARAM_DESCRIPTION
#' @param varnames PARAM_DESCRIPTION
#' @param different_ks PARAM_DESCRIPTION, Default: c(1, 3, 5, 7, 9, 11, 13)
#' @param folds PARAM_DESCRIPTION, Default: c()
#' @param do_smote PARAM_DESCRIPTION, Default: FALSE
#' @param smote_params PARAM_DESCRIPTION, Default: list(K = 5, dup_size = "balance")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}
#' @rdname makeKnn_l1o
#' @export 
#' @importFrom dplyr mutate select
makeKnn_l1o <- function(datasc, levs, varnames,
                        different_ks=c(1, 3, 5, 7, 9, 11, 13),
                        folds=c(),
                        do_smote=FALSE,
                        smote_params=list(K=5, dup_size="balance")){
  library(class)
  #library(gmodels)

  results <- list()
  datasc <- datasc %>% dplyr::mutate(class=factor(class))
  train_df_all <- datasc %>%
    dplyr::select(-class, -sample) %>%
    dplyr::select(all_of(varnames))

  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  for(k in  different_ks){
    kname = paste("K=", as.character(k), sep="", collapse="")
    results[[kname]] <- list()
    preds <- c()
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
      kires <- knn(train_df, test_df, train_labels, k = k, prob = T)
      preds  <- c(preds, kires)
      pred_probs[[i]] <- kires

    }
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



    preds2 <- knn(train_df_all, train_df_all, datasc$class, k = k, prob = T)
    results[[kname]][["preds_no_l1o"]] <- levs[preds2] %>% factor(levels=levs)
    results[[kname]][["confmat_no_l1o"]] <- confusionMatrix(results[[kname]][["preds_no_l1o"]],
                                                            factor(datasc$class),
                                                            positive = levs[2])
  }

  return(results)
}
