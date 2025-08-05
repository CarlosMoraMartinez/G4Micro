make_randomForest_l1o <- function(datasc, levs, varnames,
                                  folds=folds(),
                                  randomforest_params = randomforest_params,
                                  do_smote=FALSE,
                                  smote_params=list(K=5, dup_size="balance")
){
  library(randomForest)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))

  if(randomforest_params$balance_weights & !do_smote){
    classweights <- table(datasc$class)
    sweights <- max(classweights)/classweights[datasc$class]
  }else{
    sweights <- rep(1, nrow(datasc))
  }


  predict_tree1 <- factor()
  predict_tree1_probs <- list()


  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    train_weighs <- sweights[-i]

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
      train_weighs <- rep(1, nrow(train_df))
    }else{
      smoteData = NULL
    }
    mod_tree1 <- randomForest(x=train_df, y=train_labels, levels=levs,
                              weights = train_weighs,
                              ntree = randomforest_params$ntree,
                              mtry = randomforest_params$mtry,
                              nodesize = randomforest_params$nodesize)
    predict_tree1 <- c(predict_tree1, predict(mod_tree1, test_df))
    predict_tree1_probs[[i]] <- predict(mod_tree1, test_df, type = "prob")

  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict_tree1_probs, \(xx) xx %>% as.data.frame) %>%
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- predict_tree1_probs %>% map(as.data.frame) %>% bind_rows
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
    roc_auc <- as.numeric(roc1$auc)
  }

  mod_tree1 <- randomForest(x=df, y=datasc$class, levels=levs,
                            weights = sweights,
                            ntree = randomforest_params$ntree,
                            mtry = randomforest_params$mtry,
                            nodesize = randomforest_params$nodesize)
  predict_tree2 <- predict(mod_tree1, df)
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  return(list(confmat=confmat_tree1,
              confmat_no_l1o=confmat_tree2,
              mod=mod_tree1,
              preds=predict_tree1,
              pred_probs =probs_vector,
              preds_no_l1o=predict_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc),
              smoteData=smoteData,
              params =  randomforest_params))
}
