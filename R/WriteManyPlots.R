#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param plotlist PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 12
#' @param h PARAM_DESCRIPTION, Default: 7
#' @param separate PARAM_DESCRIPTION, Default: F
#' @param opt PARAM_DESCRIPTION, Default: list()
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname WriteManyPlots
#' @export 
WriteManyPlots <- function(plotlist, name, outdir, w=12, h=7, separate=F, opt=list()){
  if(! dir.exists(outdir)){dir.create(outdir)}
  if(separate){
    for(p in names(plotlist)){
      fname = paste0(outdir, "/", name, "_", p, ".pdf")
      pdf(fname, width=w, height=g)
      tryCatch({print(plotlist[[p]])}, error = function(x){print(getNullPlot(opt, p))})
      dev.off()
    }
  }else{
    fname = paste0(outdir, "/", name, "_all.pdf")
    pdf(fname, width=w, height=h)
    for(p in names(plotlist)){
      tryCatch({print(plotlist[[p]])}, error = function(x){print(getNullPlot(opt, p))})
    }
    dev.off()
  }
}
