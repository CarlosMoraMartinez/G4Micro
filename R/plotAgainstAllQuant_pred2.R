
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param alldf PARAM_DESCRIPTION
#' @param groupvars PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param col_var PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param names_cyt PARAM_DESCRIPTION, Default: c()
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[tidyr]{gather}}
#' @rdname plotAgainstAllQuant_pred2
#' @export 
#' @importFrom tidyr gather
plotAgainstAllQuant_pred2<- function(alldf, groupvars, name, outdir, col_var, opt, names_cyt=c()){

  plotlist <- list()
  for(grouping_var in groupvars){
    df <- alldf[, c(grouping_var, names_cyt, col_var)] %>% tidyr::gather(IL, pg_mL, names_cyt)
    oname <- paste0(outdir, name)
    plotlist[[grouping_var]] <- plotRegr(df, variable="pg_mL",
                                         x_var = grouping_var,
                                         wrap_var = "IL",
                                         col_var ="Psoriasis",
                                         fname=oname,
                                         write = F)
  }
  WriteManyPlots(plotlist, name = name, outdir =outdir, w = 12, h = 12, separate = F, opt)
  return(plotlist)
}
