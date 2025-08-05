#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param levs PARAM_DESCRIPTION
#' @param varnames PARAM_DESCRIPTION
#' @param kernel PARAM_DESCRIPTION, Default: 'linear'
#' @param SEED PARAM_DESCRIPTION, Default: 123
#' @param folds PARAM_DESCRIPTION, Default: c()
#' @param do_smote PARAM_DESCRIPTION, Default: FALSE
#' @param smote_params PARAM_DESCRIPTION, Default: list(K = 5, dup_size = 2)
#' @param balance_classes PARAM_DESCRIPTION, Default: TRUE
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
#'  \code{\link[e1071]{svm}}
#' @rdname make_svm_l1o
#' @export 
#' @importFrom dplyr select
#' @importFrom e1071 svm
make_svm_l1o <- function(datasc, levs, varnames, kernel="linear", SEED=123, folds=c(),
                         do_smote=FALSE,
                         smote_params=list(K=5, dup_size=2),
                         balance_classes=TRUE){
  library(e1071)
  datasc$class <- factor(datasc$class)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  predict1 <- factor(levels=levs)
  predict1_probs <- list()

  if(balance_classes & ! do_smote){
    class_weight <- "inverse"
  }else{
    class_weight <- rep(1, length(levs))
    names(class_weight) <- levs
  }

  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  set.seed(SEED)
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
      train_df <- smoteData %>% dplyr::select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
    }else{
      smoteData = NULL
    }

    mod <- e1071::svm(x = train_df, y = train_labels, scale=TRUE, kernel=kernel,
                      class.weights = class_weight,
                      probability = TRUE)
    predict1 <- c(predict1, predict(mod, test_df))
    predict1_probs[[i]] <- predict(mod, test_df, probability = TRUE)

  }

  confmat1 <- confusionMatrix(predict1, datasc$class, positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>%
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>%
      bind_rows #%>% pull(!!sym(levs[2]))
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
    roc_auc <- as.numeric(roc1$auc)
  }

  mod_all <- e1071::svm(x = df, y = datasc$class, scale=TRUE, kernel=kernel,
                        probability = TRUE, class.weights = class_weight)
  predict2 <- predict(mod_all, df)
  predict2_probs <- predict(mod_all, df, probability = TRUE)
  confmat2 <- confusionMatrix(predict2, datasc$class, positive = levs[2])

  mod_all_noscale <- e1071::svm(x = df[, varnames[1:2]], y = datasc$class, scale=FALSE,
                                kernel=kernel, class.weights = "inverse")
  predict2_noscale <- predict(mod_all_noscale, df[, varnames[1:2]])
  confmat2_noscale <- confusionMatrix(predict2_noscale, datasc$class, positive = levs[2])

  return(list(confmat=confmat1,
              confmat_no_l1o=confmat2,
              mod=mod_all,
              preds=predict1,
              pred_probs = probs_vector,
              pred_probs_obj = predict1_probs,
              preds_no_l1o=predict2,
              mod_noscale=mod_all_noscale,
              preds_noscale=predict2_noscale,
              confmat_noscale=confmat2_noscale,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc))
  )
}
