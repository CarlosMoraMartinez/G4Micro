#' @title Generate PCoA Plots for a Single Metadata Variable
#'
#' @description
#' This function generates PCoA (or other ordination) plots for a given phyloseq object and ordination result.
#' It creates scatter plots of the first PCoA axis (Axis 1) against additional axes specified in \code{extradims},
#' colored by a given metadata variable, and returns them combined as a single figure using \code{cowplot::plot_grid}.
#'
#' @param phobj A \code{phyloseq} object containing the microbiome data.
#' @param pcoa.bray The ordination result (e.g., output from \code{ordinate()}).
#' @param evals Eigenvalues extracted from the ordination result (typically \code{pcoa.bray$values$Eigenvalues}).
#' @param var2color Name of the metadata variable to use for coloring the plots. Default: \code{"Condition"}
#' @param name Title to be displayed on the plots, typically the distance name. Default: \code{"Bray-Curtis"}
#' @param extradims Numeric vector of additional axes (besides Axis 1) to plot. Default: \code{2:5}
#' @param labelsamples Name of the metadata column to use for labeling samples. Default: \code{"sampleID"}
#'
#' @return A combined \code{ggplot2} object (via \code{cowplot::plot_grid}) with one plot per extra axis specified in \code{extradims}.
#'
#' @details
#' This function is designed to be used internally within \code{\link{makeAllPCoAs}}.
#' For each axis in \code{extradims}, it creates a 2D plot using the first axis and the extra one.
#' Each plot includes ellipses (95% confidence), sample points, and sample labels.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[cowplot]{plot_grid}}
#'  \code{\link[phyloseq]{ordinate}},
#'  \code{\link[phyloseq]{plot_ordination}},
#' @rdname makePCoA
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom ggrepel geom_text_repel
#' @importFrom phyloseq plot_ordination
makePCoA <- function(phobj, pcoa.bray, evals,
                     var2color="Condition", name = "Bray-Curtis",
                     extradims= 2:5, labelsamples="sampleID"){
  gs <- lapply(extradims, FUN=function(axis){
    gg <- plot_ordination(phobj, pcoa.bray, color = var2color,
                          title = name, axes=c(1, axis)) +
      #coord_fixed(sqrt(evals[2] / evals[1])) +
      #scale_color_manual(values=palette2)+
      stat_ellipse(level=0.95, linetype=2, alpha = 0.8, na.rm = TRUE) +
      geom_point(size = 2) +
      ggrepel::geom_text_repel(aes_string(label = labelsamples)) +
      mytheme
    # if(!is.numeric(gg$data[, var2color])){
    #   gg <- gg + scale_color_lancet() + scale_fill_lancet()
    # }
    return(gg)
  })
  cw <- cowplot::plot_grid(plotlist=gs)
  return(cw)
}
