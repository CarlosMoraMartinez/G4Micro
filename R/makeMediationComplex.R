#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param xname PARAM_DESCRIPTION, Default: 'Edad'
#' @param yname PARAM_DESCRIPTION, Default: 'Condition_bin'
#' @param medname PARAM_DESCRIPTION, Default: 'IMC'
#' @param medname2 PARAM_DESCRIPTION, Default: 'Faecalibacterium_prausnitzii'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname makeMediationComplex
#' @export 
makeMediationComplex <- function(df, xname="Edad", yname="Condition_bin",
                                 medname="IMC", medname2="Faecalibacterium_prausnitzii"){

  eq1 = paste0(yname, " = b*", medname, " + cp*", medname2, " + fp*", xname)
  eq2 = paste0(medname, " = a*", medname2, " + d*", xname)
  eq3 = paste0(medname2, " = e*", xname)
  eqs <- paste(eq1, eq2, eq3, sep="\n", collapse="\n")

  iris.model<-specifyEquations(exog.variances=T, text=eqs)

  effects<-c('a*b', 'cp+a*b', 'd*b', 'fp+d*b', 'e*cp', 'fp+e*cp', 'd*b+fp+e*cp+e*a*b')
  nlsy.res<-bmem.sobel(df, iris.model, effects)
  return(nlsy.res)

}
