#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param var PARAM_DESCRIPTION, Default: 'Condition'
#' @param outname PARAM_DESCRIPTION, Default: 'phylumTests'
#' @param paired PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[phyloseq]{tax_glom}}, \code{\link[phyloseq]{taxa_names}}, \code{\link[phyloseq]{tax_table}}, \code{\link[phyloseq]{psmelt}}
#'  \code{\link[dplyr]{group_map}}
#'  \code{\link[broom]{reexports}}
#' @rdname getPhylumTests
#' @export 
#' @importFrom phyloseq tax_glom taxa_names tax_table psmelt
#' @importFrom dplyr group_modify
#' @importFrom broom tidy
getPhylumTests <- function(phobj, var="Condition", outname="phylumTests", paired=F){
  ## Total abundance by phylum, apparently by summing over all ASVs in phylum
  ps_phylum <- phyloseq::tax_glom(phobj, "Phylum")
  phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]

  comps <- combn(unique(unlist(sample_data(ps_phylum)[, var])), 2, simplify = F)
  signif_codes <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

  dfmelt <- phyloseq::psmelt(ps_phylum) %>%
    group_by(Phylum)

  pvals <- dfmelt %>%
    dplyr::group_modify(~ broom::tidy(wilcox.test(data=., as.formula(paste0( "Abundance ~", var)), paired = paired)))
  pvals$padj = p.adjust(pvals$p.value, method="BH")
  write_tsv(pvals, file=outname)
  return(pvals)
}
