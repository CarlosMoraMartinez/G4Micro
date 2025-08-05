#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION, Default: 'CCA'
#' @param outdir PARAM_DESCRIPTION, Default: 'CCA'
#' @param variableList PARAM_DESCRIPTION, Default: list(c("Condition"))
#' @param dist_type PARAM_DESCRIPTION, Default: 'bray'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname makeConstrainedOrdinationAll
#' @export 
makeConstrainedOrdinationAll <- function(phobj, opt, name="CCA", outdir= "CCA", variableList = list(c("Condition")),  dist_type="bray"){
  plots <- list()
  i = 1
  for(vars in variableList){
    if(length(vars) ==1 ){
      plots[[i]] <- makeConstrainedOrdinationSingleVar(phobj, variables =vars,  dist_type=dist_type)
    }else{
      plots[[i]] <- makeConstrainedOrdination(phobj, variables =vars,  dist_type=dist_type)
    }
    i <- i+1
  }
  names(plots)<- sapply(variableList, paste, sep="-", collapse="-")
  WriteManyPlots(plots, name, outdir, opt=opt)
  return(plots)
}
