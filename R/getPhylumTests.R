#' @title Perform Wilcoxon Tests on Phylum Abundances
#' @description Conducts pairwise Wilcoxon rank-sum tests on phylum-level abundance data extracted from a phyloseq object. Returns and saves a table with p-values and adjusted p-values (FDR).
#' @param phobj A \code{phyloseq} object containing microbiome data.
#' @param var A character string indicating the sample variable to group by for comparisons. Default: 'Condition'
#' @param outname Output filename for the results table. Default: 'phylumTests'
#' @param paired Logical. Should the test be paired? Default: FALSE
#' @return A data frame with Wilcoxon test results for each phylum, including adjusted p-values.
#' @details This function performs non-parametric testing for differential abundance at the phylum level. It uses Wilcoxon rank-sum tests (or signed-rank if \code{paired=TRUE}) and adjusts p-values using the Benjamini-Hochberg method.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[phyloseq]{tax_glom}},
#'  \code{\link[phyloseq]{taxa_names}},
#'  \code{\link[phyloseq]{tax_table}},
#'  \code{\link[phyloseq]{psmelt}}
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
