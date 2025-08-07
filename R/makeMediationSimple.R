#' @title Run a Mediation Analysis
#' @description
#' Performs a basic mediation analysis using the `bmem` and `sem` packages. The function estimates
#' the direct, indirect, and total effects of an independent variable on a dependent variable through a mediator.
#'
#' @param df Data frame. The dataset containing the variables for mediation analysis.
#' @param xname Character. Name of the independent (predictor) variable in \code{df}.
#' @param yname Character. Name of the dependent (outcome) variable in \code{df}.
#' @param medname Character. Name of the mediator variable in \code{df}.
#'
#' @return An object of class \code{bmem}, containing the results of the mediation analysis,
#' including estimated paths, indirect and total effects, and significance.
#'
#' @details
#' The function fits a two-equation mediation model:
#' \itemize{
#'   \item \code{medname = a * xname} (effect of the independent variable on the mediator)
#'   \item \code{yname = b * medname + cp * xname} (effect of the mediator and the direct effect of the independent variable on the outcome)
#' }
#' It then uses the Sobel test to estimate the indirect effect (a*b) and total effect (cp + a*b).
#' The model is specified using `specifyEquations` from the `sem` package and estimated with `bmem.sobel` from the `bmem` package.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname makeMediationSimple
#' @export
#' @importFrom bmem bmem.sobel
#' @importFrom sem specifyEquations
makeMediationSimple <- function(df, xname, yname, medname){
  library(bmem)
  library(sem)

  eq1 = paste0(yname, " = b*", medname, " + cp*", xname)
  eq2 = paste0(medname, " = a*", xname)
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")

  iris.model<-specifyEquations(exog.variances=T, text=eqs)

  effects<-c('a*b', 'cp+a*b')
  nlsy.res<-bmem.sobel(df, iris.model, effects)
  #nlsy.res<-bmem.bs(df, iris.model, effects)
  return(nlsy.res)

}
