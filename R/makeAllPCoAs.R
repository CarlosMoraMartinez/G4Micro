#' @title Create and Save Multiple PCoA Plots from a Phyloseq Object
#' @description
#' This function computes Principal Coordinates Analysis (PCoA) or other ordinations on a given phyloseq object,
#' generates plots colored by specified metadata variables, and optionally saves them as multi-panel PDFs.
#'
#' @param phobj A \code{phyloseq} object containing the microbiome data.
#' @param outdir Output directory where the plots will be saved.
#' @param method Ordination method to use (e.g., \code{"PCoA"}, \code{"NMDS"}). Default: \code{"PCoA"}
#' @param name A prefix name for saved files and plots. Default: \code{"PCoAs"}
#' @param dist_type Dissimilarity measure to use for ordination (e.g., \code{"bray"}, \code{"jaccard"}). Default: \code{"bray"}
#' @param dist_name Pretty name of the distance metric to display in plot titles. Default: \code{"Bray-Curtis"}
#' @param vars2plot A character vector of metadata variables to use for coloring the plots. If empty, all available variables (excluding "Codigo" and "Num_paciente") will be used. Default: \code{c()}
#' @param extradims Numeric vector of additional PCoA axes to plot beyond the first two (e.g., 2:5). Default: \code{2:5}
#' @param create_pdfs Logical; if \code{TRUE}, plots are saved to multi-page PDF files. If \code{FALSE}, only a single plot is returned. Default: \code{MULTI_PAGE_PDFS}
#' @param labelsamples Name of the metadata column to use for sample labels in the plots. Default: \code{"sampleID"}
#' @param w Width of the saved plots (in inches). Default: \code{12}
#' @param h Height of the saved plots (in inches). Default: \code{5}
#'
#' @return A named list of ggplot2 plot objects for each metadata variable. If \code{create_pdfs = FALSE}, only one plot is returned in the list under the name "Condition".
#'
#' @details
#' This function performs ordination (default: PCoA) on the provided phyloseq object using the specified distance metric.
#' It creates scatter plots of the first two ordination axes (and optionally additional dimensions),
#' colored by sample metadata variables. When \code{create_pdfs = TRUE}, all plots are saved as a multi-panel PDF in the given output directory.
#' If no variables are specified in \code{vars2plot}, it will automatically use all metadata variables except \code{"Codigo"} and \code{"Num_paciente"}.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  pcoas <- makeAllPCoAs(phobj,
#'                         outdir = "results/pcoas",
#'                         vars2plot = c("Condition", "Sex"),
#'                         dist_type = "bray")
#'   pcoas$Condition  # Access the plot for the variable "Condition"
#'  }
#' }
#' @seealso
#'  \code{\link[phyloseq]{ordinate}},
#'  \code{\link[phyloseq]{plot_ordination}},
#'  \code{\link{makePCoA}}
#'  \code{\link{WriteManyPlots}}
#'
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
