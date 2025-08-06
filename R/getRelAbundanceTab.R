#' @title Get Relative Abundance and Prevalence Table
#' @description Calculates prevalence and total abundance of taxa from a phyloseq object and returns a data frame with taxonomy information.
#' @param phobj A phyloseq object containing microbiome count data and taxonomy.
#' @return A data frame with columns:
#' \itemize{
#'   \item Prevalence: number of samples where taxa is present (count > 0)
#'   \item TotalAbundance: sum of counts across all samples
#'   \item Taxonomy columns from the phyloseq object's tax_table
#'   \item relative_prevalence: prevalence divided by total number of samples
#' }
#' @details
#' Prevalence is calculated as the count of samples in which each taxon has a non-zero abundance.
#' Total abundance is the sum of counts across all samples.
#' Relative prevalence is the proportion of samples where the taxon is present.
#' This function works for phyloseq objects with taxa as rows or columns.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[phyloseq]{otu_table}},
#'  \code{\link[phyloseq]{taxa_sums}}
#' @rdname getRelAbundanceTab
#' @export
#' @importFrom phyloseq otu_table taxa_sums tax_table nsamples taxa_are_rows
getRelAbundanceTab <- function(phobj){
  ottmp <- phyloseq::otu_table(phobj)

  pre_prevalence <- apply(X = ottmp,
                          MARGIN = ifelse(taxa_are_rows(phobj), yes = 1, no = 2),
                          FUN = function(x){sum(x > 0)})
  pre_prevalence = data.frame(Prevalence = pre_prevalence,
                              TotalAbundance = phyloseq::taxa_sums(phobj),
                              tax_table(phobj),
                              relative_prevalence = pre_prevalence/ nsamples(phobj)
  )
  return(pre_prevalence)
}
