make_catboost_l1o <- function(datasc, levs, varnames,
                              folds=folds(),
                              catboost_params = catboost_params,
                              do_smote=FALSE,
                              smote_params=list(K=5, dup_size="balance")
){
  library(catboost)
  datasc$class <- factor(datasc$class, levels=levs)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  predict_probs = numeric(0)

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
      train_pool <- catboost.load_pool(data = train_df, label = as.integer(train_labels == levs[2]))
    }else{
      train_pool <- catboost.load_pool(data = train_df, label = as.integer(train_labels)-1)
    }
    test_pool <- catboost.load_pool(data = test_df)

    if(catboost_params$bootstrap_type == "Bernoulli"){
      model <- catboost.train(learn_pool = train_pool, params = list(
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
      model <- catboost.train(learn_pool = train_pool, params = list(
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
    pred_prob <- catboost.predict(model, test_pool, prediction_type = "Probability")

    if(length(levs) == 2){
      predict_probs <- c(predict_probs, pred_prob)
    }else{
      predict_probs <- rbind(predict_probs, pred_prob)
    }

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

  train_pool <- catboost.load_pool(data = df, label = as.integer(datasc$class == levs[2]))
  test_pool <- catboost.load_pool(data = df)

  if(catboost_params$bootstrap_type == "Bernoulli"){
    model2 <- catboost.train(learn_pool = train_pool, params = list(
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
    model2 <- catboost.train(learn_pool = train_pool, params = list(
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
  predict_tree2 <- catboost.predict(model2, test_pool, prediction_type = "Probability")
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
