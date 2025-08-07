#' @title Perform Mediation Analysis with Multiple Predictors Merged into a Single Model
#' @description Runs a simple mediation model where multiple independent variables (`xnames`)
#' affect a single mediator (`medname`), which in turn affects the dependent variable (`yname`).
#' The function constructs the equations dynamically and calculates direct, indirect, and total effects using the Sobel test.
#' @param df A data frame containing the variables involved in the mediation model.
#' @param xnames A character vector with the names of the independent variables.
#' @param yname A character string indicating the name of the dependent (outcome) variable.
#' @param medname A character string indicating the name of the mediator variable.
#' @return A data frame with mediation effect estimates, standard errors, z-scores, p-values, and matched x variable labels for each path.
#' @details
#' This function supports multiple independent variables by dynamically constructing the mediation equations. It uses the `bmem.sobel` function from the `bmem` package to estimate mediation effects and returns both the indirect and total effects per variable. Useful for simple mediation power calculations or exploratory analyses.
#'
#' The function adds a `x_labels` column to indicate which independent variable each effect corresponds to.
#'
#' Note: This function is designed for continuous or dichotomous outcomes and assumes linear relationships.
#'@examples
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
