makeKnn <- function(datasc, levs, nvars, different_ks=c(1, 3, 5, 7, 9, 11, 13)){
  library(class)
  #library(gmodels)

  test_pred <- list()
  conf_matrices_knn <- list()
  train_df <- datasc %>% dplyr::select(-class, -sample) %>% dplyr::select(all_of(varnames))
  train_labels <- datasc$class %>% factor(levels=levs)

  for(k in  different_ks){
    kname = paste("K=", as.character(k), sep="", collapse="")
    test_pred[[kname]] <- knn(train_df, train_df, train_labels, k = k, prob = T)
    conf_matrices_knn[[kname]] <- confusionMatrix(test_pred[[kname]],
                                                  train_labels,
                                                  positive = levs[2])
  }

  return(list(confmats=conf_matrices_knn, mods=test_pred))
}
