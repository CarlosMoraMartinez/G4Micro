makeMediationSimple_mergedMediat <- function(df, x_name, y_name, mednames){

  a_params <- paste("a", 1:length(mednames), sep="")
  b_params <- paste("b", 1:length(mednames), sep="")

  eq1 = paste(y_name, " = ", b_params, "*", mednames, " + ", "cp*", x_name, sep="") %>%
    paste(collapse = "\n")
  eq2 = paste0(mednames, " = ", a_params, "*", x_name) %>%
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")

  iris.model<-specifyEquations(exog.variances=T, text=eqs)

  effects<-c(paste(a_params, "*", b_params, sep="" ), paste('cp+', a_params,'*', b_params, sep=""))
  nlsy.res<-bmem.sobel(df, iris.model, effects)

  oldnames2 <- rownames(nlsy.res$estimates) %>%
    sapply(\(x){y<-strsplit(x, "\\*|\\+",perl = TRUE)[[1]];y[length(y)]})
  # newnames <- ifelse(oldnames %in% a_params,
  #        paste(oldnames, xnames[match(oldnames, a_params)], sep="_"),
  #        ifelse(oldnames %in% cp_params,
  #               paste(oldnames, xnames[match(oldnames, cp_params)], sep="_"),
  #               oldnames))
  resdf <- nlsy.res$estimate %>%
    rownames_to_column("param")
  resdf$x_labels <- ifelse(oldnames2 %in% a_params,
                           mednames[match(oldnames2, a_params)],
                           ifelse(oldnames2 %in% b_params,
                                  mednames[match(oldnames2, b_params)],
                                  gsub("V\\[|\\]", "", oldnames2)))

  #rownames(nlsy.res$estimates) <- newnames
  return(resdf)
}
