#' @title Plot Prevalence vs Abundance by Phylum
#' @description Creates a scatter plot of total abundance vs relative prevalence for taxa grouped by Phylum using data from a phyloseq object.
#' @param phobj A phyloseq object containing microbiome data.
#' @param outname File name for the output PDF plot. Default: 'phylumBarplot.pdf'
#' @param height Height of the output PDF plot in inches. Default: 10
#' @param width Width of the output PDF plot in inches. Default: 12
#' @return A `ggplot` object showing prevalence vs abundance faceted by Phylum.
#' @details The plot displays taxa's total abundance (log scale) against their prevalence across samples, with a horizontal reference line at 0.05 prevalence.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link{getRelAbundanceTab}}
#' @rdname plotPrevalenceVsAbundance
#' @export
#' @importFrom ggpubr theme_pubclean
plotPrevalenceVsAbundance <- function(phobj, outname="phylumBarplot.pdf", height=10, width=12){
  pre_prevalence <- getRelAbundanceTab(phobj)
  g1 <- ggplot(pre_prevalence, aes(TotalAbundance, relative_prevalence, color = Phylum)) +
    # Agregamos una línea para nuestro umbral
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") +
    ylab("Prevalence/num. samples") +
    facet_wrap(~Phylum) +
    theme_pubclean() +
    guides(color = FALSE)
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}
