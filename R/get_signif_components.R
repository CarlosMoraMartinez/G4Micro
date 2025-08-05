get_signif_components <- function(datasc, levs){
  df <- datasc %>% dplyr::select(-sample)
  varnames <- names(df)[names(df)!="class"]
  res <- data.frame()
  ps <- c()
  for(i in varnames){
    df$aux <- df[, i]
    mod_glm <- glm(class ~ aux, data=df, family = binomial)
    ps <- c(ps, summary(mod_glm)$coefficients[2, 4])
    predict_glm1 <- predict(mod_glm, df)
    predict_glm1 <- ifelse(predict_glm1 > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
    confmat <- confusionMatrix(predict_glm1, datasc$class, positive = levs[2])
    auxdf <- data.frame(var=i,
                        pval=summary(mod_glm)$coefficients[2, 4],
                        Accuracy=confmat$overall["Accuracy"],
                        Sensitivity = confmat$byClass["Sensitivity"],
                        Specificity = confmat$byClass["Specificity"],
                        PPV = confmat$byClass["Pos Pred Value"],
                        NPV = confmat$byClass["Neg Pred Value"]
    )
    res <- rbind(res, auxdf)
  }
  return(res %>% dplyr::arrange(pval))
}
