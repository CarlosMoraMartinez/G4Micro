#' @title Plot Boxplots of Phylum-Level Abundances
#' @description
#' Plots boxplots of phylum-level abundances using a phyloseq object, faceted by phylum,
#' with optional paired comparisons between groups. Significance is assessed using Wilcoxon tests.
#' @param phobj A phyloseq object containing microbiome data.
#' @param var A sample metadata variable used to group the boxplots, Default: 'Condition'
#' @param outname Name of the output PDF file, Default: 'phylumBarplot.pdf'
#' @param height Height (in inches) of the output PDF, Default: 10
#' @param width Width (in inches) of the output PDF, Default: 8
#' @param paired Logical indicating if paired Wilcoxon tests should be used, Default: FALSE
#' @return A ggplot2 object representing the boxplots.
#' @details
#' This function summarizes abundances at the phylum level, then creates boxplots grouped by the
#' specified metadata variable. Significance comparisons between groups are performed using
#' `stat_compare_means` from the ggpubr package.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plotPhylumBoxplots(my_phyloseq_object)
#'  }
#' }
#' @seealso
#'  \code{\link[phyloseq]{tax_glom}},
#'  \code{\link[phyloseq]{taxa_names}},
#'  \code{\link[phyloseq]{tax_table}},
#'  \code{\link[phyloseq]{psmelt}}
#'  \code{\link[ggpubr]{stat_compare_means}},
#'  \code{\link[ggplot2]{ggplot}},
#'  \code{\link[ggplot2]{geom_boxplot}},
#'  \code{\link[ggplot2]{geom_jitter}},
#'  \code{\link[ggplot2]{geom_line}},
#'  \code{\link[ggplot2]{facet_wrap}}
#' @rdname plotPhylumBoxplots
#' @export
#' @importFrom phyloseq tax_glom taxa_names tax_table psmelt
#' @importFrom ggpubr stat_compare_means
#' @importFrom ggpubr theme_pubclean
plotPhylumBoxplots <- function(phobj, var="Condition", outname="phylumBarplot.pdf", height=10, width=8, paired=F){
  ## Total abundance by phylum, apparently by summing over all ASVs in phylum
  ps_phylum <- phyloseq::tax_glom(phobj, "Phylum")
  phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]

  sample_data(ps_phylum)[, var] <- unlist(sample_data(ps_phylum)[, var]) %>% as.character
  comps <- combn(unique(unlist(sample_data(ps_phylum)[, var])), 2, simplify = F)
  signif_codes <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

  g1 <- phyloseq::psmelt(ps_phylum) %>%
    ggplot(data = ., aes_string(x = var, y = "Abundance")) +
    facet_wrap(~ OTU, scales = "free")+
    geom_boxplot(outlier.shape  = NA) +
    geom_jitter(aes(color = OTU), height = 0, width = .2) +
    ggpubr::stat_compare_means(method="wilcox.test", comparisons = comps,
                       symnum.args = signif_codes, paired=paired) +
    labs(x = var, y = "Abundance\n") +
    ggpubr::theme_pubclean() +
    guides(color = FALSE)
  if(paired){
    g1 <- g1 + geom_line(aes(group = pacienteID), color = "gray20", linetype=1, size=0.1, alpha=0.25)
  }
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}
