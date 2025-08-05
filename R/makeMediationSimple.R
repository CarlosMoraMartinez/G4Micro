#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param xname PARAM_DESCRIPTION
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
#' @rdname makeMediationSimple
#' @export 
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
