#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param levs PARAM_DESCRIPTION
#' @param varnames PARAM_DESCRIPTION
#' @param folds PARAM_DESCRIPTION, Default: folds()
#' @param xgboost_params PARAM_DESCRIPTION, Default: xgboost_params
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
#' @rdname make_xgboost_l1o
#' @export 
#' @importFrom dplyr select
make_xgboost_l1o <- function(datasc, levs, varnames,
                             folds=folds(),
                             xgboost_params = xgboost_params,
                             do_smote=FALSE,
                             smote_params=list(K=5, dup_size="balance")
){
  library(xgboost)
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
