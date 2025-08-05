#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param levs PARAM_DESCRIPTION
#' @param varnames PARAM_DESCRIPTION
#' @param folds PARAM_DESCRIPTION, Default: c()
#' @param do_smote PARAM_DESCRIPTION, Default: FALSE
#' @param smote_params PARAM_DESCRIPTION, Default: list(K = 5, dup_size = 2)
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}
#'  \code{\link[assertthat]{assert_that}}
#' @rdname make_glm_l1o_multiclass
#' @export 
#' @importFrom dplyr select mutate
#' @importFrom assertthat assert_that
make_glm_l1o_multiclass <- function(datasc, levs, varnames, folds=c(),
                                    do_smote=FALSE,
                                    smote_params=list(K=5, dup_size=2)){
  library(nnet)
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

    mod_glm <- multinom(formula, data=train_df)
    newpred <- predict(mod_glm, test_df, type="prob")
    predict_glm1 <- rbind(predict_glm1, newpred)

  }

  predict1<- colnames(predict_glm1)[apply(predict_glm1, MAR=1, which.max)] %>% factor
  confmat1 <- confusionMatrix(predict1, factor(datasc$class))

  assertthat::assert_that(all(colnames(predict_glm1) == levels(datasc$class)))
  roc1 <- multiclass.roc(response=datasc$class, predictor=predict_glm1)
  roc_auc <- as.numeric(roc1$auc)

  mod_all <- multinom(formula, data=datasc, family = binomial)
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
