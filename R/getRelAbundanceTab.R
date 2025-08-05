#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[phyloseq]{otu_table}}, \code{\link[phyloseq]{taxa_sums}}
#' @rdname getRelAbundanceTab
#' @export 
#' @importFrom phyloseq otu_table taxa_sums
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
