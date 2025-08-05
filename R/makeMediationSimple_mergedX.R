#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param xnames PARAM_DESCRIPTION
#' @param yname PARAM_DESCRIPTION
#' @param medname PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname makeMediationSimple_mergedX
#' @export 
makeMediationSimple_mergedX <- function(df, xnames, yname, medname){

  a_params <- paste("a", 1:length(xnames), sep="")
  cp_params <- paste("cp", 1:length(xnames), sep="")

  eq1 = paste(yname, " = b*", medname, " + ", cp_params, "*", xnames, sep="") %>%
    paste(collapse = "\n")
  eq2 = paste0(medname, " = ", a_params, "*", xnames) %>%
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")

  iris.model<-specifyEquations(exog.variances=T, text=eqs)

  effects<-c(paste(a_params, "*b", sep="" ), paste(cp_params, '+', a_params,'*b', sep=""))
  nlsy.res<-bmem.sobel(df, iris.model, effects)
  oldnames2 <- rownames(nlsy.res$estimates) %>%
    sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1])
  # newnames <- ifelse(oldnames %in% a_params,
  #        paste(oldnames, xnames[match(oldnames, a_params)], sep="_"),
  #        ifelse(oldnames %in% cp_params,
  #               paste(oldnames, xnames[match(oldnames, cp_params)], sep="_"),
  #               oldnames))
  resdf <- nlsy.res$estimate %>%
    rownames_to_column("param")
  resdf$x_labels <- ifelse(oldnames2 %in% a_params,
                           xnames[match(oldnames2, a_params)],
                           ifelse(oldnames2 %in% cp_params,
                                  xnames[match(oldnames2, cp_params)],
                                  gsub("V\\[|\\]", "", oldnames2)))
  #rownames(nlsy.res$estimates) <- newnames
  return(resdf)
}
