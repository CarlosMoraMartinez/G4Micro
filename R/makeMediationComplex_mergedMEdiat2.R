#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param xname PARAM_DESCRIPTION, Default: 'Edad'
#' @param yname PARAM_DESCRIPTION, Default: 'Condition_bin'
#' @param medname PARAM_DESCRIPTION, Default: 'IMC'
#' @param mednames2 PARAM_DESCRIPTION, Default: c("Faecalibacterium_prausnitzii")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stringr]{str_starts}}
#' @rdname makeMediationComplex_mergedMEdiat2
#' @export 
#' @importFrom stringr str_starts
makeMediationComplex_mergedMEdiat2 <- function(df, xname="Edad",
                                               yname="Condition_bin",
                                               medname="IMC",
                                               mednames2=c("Faecalibacterium_prausnitzii")){
  C_CTRL = C_CTRL2

  a_params <- paste("a", 1:length(mednames2), sep="")
  cp_params <- paste("cp", 1:length(mednames2), sep="")
  e_params <- paste("e", 1:length(mednames2), sep="")

  eq1 = paste(yname, " = b*", medname, " + ", cp_params, "*", mednames2, " + fp*", xname, sep="") %>%
    paste(collapse = "\n")
  eq2 = paste0(medname, " = ", a_params, "*", mednames2, " + d*", xname) %>%
    paste(collapse = "\n")
  eq3 = paste0(mednames2, " = ", e_params, "*", xname)%>%
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, eq3, sep="\n", collapse="\n")

  iris.model<-specifyEquations(exog.variances=T, text=eqs)

  effects<-c(paste(a_params, "*b", sep="" ),
             paste(cp_params, '+', a_params,'*b', sep=""),
             'fp+d*b',
             paste(e_params, '*', cp_params, sep=""),
             paste('fp + ', e_params, '*', cp_params, sep=""),
             paste('d*b + fp + ', e_params, '*', cp_params, ' + ', e_params, '*',a_params, "*b", sep="")
  )
  nlsy.res<-bmem.sobel(df, iris.model, effects)

  resdf <- nlsy.res$estimate %>%
    rownames_to_column("param")
  resdf$x_labels <- resdf$param %>%
    sapply(\(x){
      y <-strsplit(x, "\\*|\\+",perl = TRUE)[[1]] %>% gsub(" ", "", .) %>% sort

      if(any(stringr::str_starts(y, "a"))){
        pthis <- y[stringr::str_starts(y, "a")]
        return(paste(mednames2[a_params==pthis], collapse=":"))
      }else if(any(stringr::str_starts(y, "cp"))){
        pthis <- y[stringr::str_starts(y, "cp")]
        return(paste(mednames2[cp_params==pthis], collapse=":"))
      }else if(any(stringr::str_starts(y, "e"))){
        pthis <- y[stringr::str_starts(y, "e")]
        return(paste(mednames2[e_params==pthis], collapse=":"))
      }else{
        return(gsub("V\\[|\\]", "", x))
      }

    })
  return(resdf)
}
