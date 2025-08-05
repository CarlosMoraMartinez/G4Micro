#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'PCoA'
#' @param name PARAM_DESCRIPTION, Default: 'PCoAs'
#' @param dist_type PARAM_DESCRIPTION, Default: 'bray'
#' @param dist_name PARAM_DESCRIPTION, Default: 'Bray-Curtis'
#' @param vars2plot PARAM_DESCRIPTION, Default: c()
#' @param extradims PARAM_DESCRIPTION, Default: 2:5
#' @param create_pdfs PARAM_DESCRIPTION, Default: MULTI_PAGE_PDFS
#' @param labelsamples PARAM_DESCRIPTION, Default: 'sampleID'
#' @param w PARAM_DESCRIPTION, Default: 12
#' @param h PARAM_DESCRIPTION, Default: 5
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname makeAllPCoAs
#' @export 
makeAllPCoAs <- function(phobj, outdir,
                         method = "PCoA",
                         name="PCoAs",
                         dist_type = "bray",
                         dist_name = "Bray-Curtis",
                         vars2plot = c(),
                         extradims = 2:5,
                         create_pdfs = MULTI_PAGE_PDFS,
                         labelsamples="sampleID", w=12, h=5){
  #palette2 <- RColorBrewer::brewer.pal(n = 9, name = 'Set1')[8:9]
  pcoa.bray <- ordinate(phobj, method = method, distance = dist_type)

  evals <- pcoa.bray$values$Eigenvalues

  if(create_pdfs){
    # Plot all the variables and save to PDF
    if(length(vars2plot) == 0){
      vars2plot <- phobj %>% sample_data() %>% names()
      vars2plot <- vars2plot[! vars2plot %in% c("Codigo", "Num_paciente")]
    }
    all_pcoas_plots <- lapply(vars2plot, FUN=function(vv, phobj, pcoa.bray, evals){
      makePCoA(phobj, pcoa.bray, evals, vv, dist_name, extradims, labelsamples = labelsamples)
    },phobj, pcoa.bray, evals)
    names(all_pcoas_plots) <- vars2plot

    WriteManyPlots(all_pcoas_plots, name, outdir, w=w, h=h, separate=F, opt)

  }else{
    #Plot only one variable and return the plot without saving it
    all_pcoas_plots <- list("Condition"=makePCoA(phobj, pcoa.bray, evals, vars2plot[1], dist_name, extradims))
  }

  return(all_pcoas_plots)
}
