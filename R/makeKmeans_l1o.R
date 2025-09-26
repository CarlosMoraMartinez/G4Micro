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
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{summarise}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{mutate_all}}
#' @rdname makeKmeans_l1o
#' @export
#' @importFrom dplyr select mutate summarise arrange mutate_all
#' @importFrom UBL SmoteClassif
#' @importFrom pROC roc multiclass.roc
#' @importFrom caret confusionMatrix
#' @importFrom stats as.formula
makeKmeans_l1o <- function(datasc, levs, varnames, SEED=123, folds=c(),
                           do_smote=FALSE,
                           smote_params=list(K=5, dup_size="balance")){
  set.seed(SEED)
  train_df_all <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))

  if(length(folds)==0){
    folds <- 1:nrow(datasc)
    reorder_samples <- FALSE
  }else{
    reorder_samples <- TRUE
  }

  preds <- c()
  all_dists <- list()
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

    mod_kmeans <- kmeans(train_df, centers=length(levs), iter.max = 100, nstart=100)

    assign_class <- table(mod_kmeans$cluster, train_labels)
    cents <- mod_kmeans$centers
    rownames(cents) <- paste0("C_", rownames(cents))
    train_dists <- dist(rbind(cents, train_df %>% as.matrix )) %>%
      as.matrix %>%
      as.data.frame %>%
      rownames_to_column("sample") %>%
      dplyr::select("sample", all_of(rownames(cents))) %>%
      filter(! sample %in% rownames(cents)) %>%
      dplyr::mutate(class = train_labels) %>%
      group_by(class) %>%
      dplyr::summarise(across(all_of(rownames(cents)), .fns = mean)) %>%
      rowwise() %>%
      dplyr::mutate(min_dist = min(c_across(all_of(rownames(cents))))) %>%
      ungroup() %>%
      dplyr::arrange(min_dist)

    free_groups <- rownames(cents)
    for(ir in 1:nrow(train_dists)){
      newgr <- free_groups[which.min(train_dists[ir, free_groups] %>% as.vector %>% unlist)]
      train_dists$cluster[ir] <- newgr
      free_groups <- free_groups[free_groups != newgr]
    }

    rownames(cents) <- train_dists$class[match(rownames(cents), train_dists$cluster)]

    dists <- dist(rbind(cents, test_df)) %>% as.matrix
    assign_class <- rownames(cents)

    if(length(i)> 1){ # already there
      pred_i <- apply(dists[rownames(test_df),  assign_class], MAR=1, \(x) assign_class[which.min(x)])
      for(ii in i) all_dists[[ii]] <- dists[as.character(ii),  assign_class]
    }else{
      pred_i <- assign_class[which.min(dists[rownames(test_df),  assign_class])]
      all_dists[[i]] <- dists[rownames(test_df),  assign_class]
    }
    preds  <- c(preds, pred_i)
  }
  if(reorder_samples){
    preds <- preds[as.character(1:nrow(datasc))] ## check
  }
  preds <- factor(preds, levels=levels(datasc$class))
  confmat_kmeans <- confusionMatrix(preds, datasc$class, positive = levs[2])
  if(length(levs) == 2){
    all_dists <- all_dists %>% bind_rows
    probs <- 1 - (all_dists %>% pull(!!sym(levs[2])))/rowSums(all_dists)
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs)
  }else{
    probs <- all_dists %>% bind_rows %>%
      dplyr::mutate_all(\(x) 1 - x/rowSums(.))
    roc1 <- multiclass.roc(response=datasc$class, predictor=as.matrix(probs))
  }

  mod_kmeans_all <- kmeans(train_df_all, centers=length(levs), iter.max = 100, nstart=100)
  predict_kmeans_nol1o <-levels(datasc$class)[mod_kmeans_all$cluster] %>% factor(levels=levs)
  confmat_kmeans_nol1o <- confusionMatrix(predict_kmeans_nol1o, datasc$class, positive = levs[2])
  if(length(levs) == 2 & confmat_kmeans_nol1o$overall["Accuracy"] < 0.5){
    predict_kmeans_nol1o <-rev(levels(datasc$class))[mod_kmeans_all$cluster] %>% factor(levels=levs)
    confmat_kmeans_nol1o <- confusionMatrix(predict_kmeans_nol1o, datasc$class, positive = levs[2])
  }
  return(list(confmat = confmat_kmeans,
              confmat_no_l1o=confmat_kmeans_nol1o,
              preds=preds,
              pred_probs=probs,
              preds_no_l1o=predict_kmeans_nol1o,
              mod=mod_kmeans_all,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc))
  )
}
