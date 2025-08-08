#' @export
NFOLDS_DEFAULT = 0

#' @export
randomforest_params_default = list(ntree = 500,
                           mtry = 1,
                           nodesize = 1,
                           balance_weights = TRUE)

#' @export
xgboost_params_default =  list(learning_rate=0.3,
                       max_depth=2,
                       nrounds =50,
                       min_child_weight=1,
                       subsample =0.6,
                       colsample_bytree =1,
                       reg_lambda =3,
                       reg_alpha =0,
                       gamma = 0,
                       nthread=1,
                       objective= "binary:logistic",
                       balance_weights = TRUE
)

#' @export
catboost_params_default <- list(
  iterations = 100,
  learning_rate = 0.05,
  depth = 2,
  loss_function = "Logloss",
  eval_metric = "AUC",
  random_seed = 123,
  use_best_model = TRUE,
  od_type = "Iter",
  od_wait = 20,
  verbose = FALSE,
  thread_count = 1,
  balance_weights = TRUE,
  bootstrap_type = "Bayesian",
  l2_leaf_reg = 3,
  subsample = 0.6,  # only used if bootstrap type ="Bernouilli"
  grow_policy = "Depthwise",
  auto_class_weights = "Balanced"
)

#' @export
smote_params_default = list(K=5, dup_size="balance")

#' @export
randomforest_params_mult_default = list(ntree = 500,
                                mtry = 4,
                                nodesize = 5,
                                balance_weights = TRUE)

#' @export
xgboost_params_mult_default =  list(learning_rate=0.3,
                            max_depth=2,
                            nrounds =30,
                            min_child_weight=1,
                            subsample =1,
                            colsample_bytree =0.6,
                            gamma = 0,
                            reg_lambda =1,
                            reg_alpha =0,
                            nthread=1,
                            objective= "multi:softprob",
                            num_class = 4,
                            balance_weights = TRUE
)

#' @export
catboost_params_mult_default <- list(
  iterations = 100,
  learning_rate = 0.05,
  depth = 4,
  loss_function = "MultiClass",
  eval_metric = "MultiClass",
  random_seed = 123,
  use_best_model = TRUE,
  od_type = "Iter",
  od_wait = 20,
  verbose = FALSE,
  thread_count = 4,
  balance_weights = TRUE, #not used anymore
  bootstrap_type = "Bernoulli",
  l2_leaf_reg = 3,
  subsample = 0.6,  # only if bootstrap type ="Bernouilli"
  grow_policy = "Depthwise",
  auto_class_weights = "Balanced"
)
