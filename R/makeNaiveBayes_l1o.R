#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param levs PARAM_DESCRIPTION
#' @param varnames PARAM_DESCRIPTION
#' @param SEED PARAM_DESCRIPTION, Default: 123
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
#'  \code{\link[dplyr]{select}}
#' @rdname makeNaiveBayes_l1o
#' @export 
#' @importFrom dplyr select
makeNaiveBayes_l1o <- function(datasc, levs, varnames,
                               SEED=123, folds=c(),
                               do_smote=FALSE,
                               smote_params=list(K=5, dup_size="balance")){
  library(e1071)
  set.seed(SEED)

  predict_bayes1 <- factor()
  predict_bayes1_probs <- list()
  df <- datasc %>% dplyr::select(-class, -sample) %>% dplyr::select(all_of(varnames))
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
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
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
    }else{
      smoteData = NULL
    }

    mod_bayes2 <- naiveBayes(train_df, train_labels, laplace = 0)
    predict_bayes1 <- c(predict_bayes1, predict(mod_bayes2, test_df))
    predict_bayes1_probs[[i]] <- predict(mod_bayes2, test_df, type="raw")

  }
  confusionMatrix_bayes1 <- confusionMatrix(predict_bayes1, datasc$class,
                                            positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict_bayes1_probs, \(xx) xx %>% as.data.frame) %>%
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- map(predict_bayes1_probs, \(xx) xx %>% as.data.frame) %>%
      bind_rows
    #probs_vector2 <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>%
    #   map2(datasc$class, \(x, nn) x[as.character(nn)]) %>% unlist
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
  }
  modwithall <- naiveBayes(df, datasc$class, laplace = 0)
  predict_bayes2 <- predict(modwithall, df)
  confMatrix_bayes2_nol1o <- confusionMatrix(predict_bayes2, datasc$class,
                                             positive = levs[2])
  return(list(confmat=confusionMatrix_bayes1,
              confmat_no_l1o=confMatrix_bayes2_nol1o,
              preds=predict_bayes1,
              pred_probs=probs_vector,
              preds_no_l1o=predict_bayes2,
              mod=modwithall,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc))
  )
}
