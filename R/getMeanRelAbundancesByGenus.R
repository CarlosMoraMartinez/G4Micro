#' @title Calculate Mean Prevalence of ASVs by Genus
#' @description Computes the mean and total prevalence of Amplicon Sequence Variants (ASVs) grouped by Genus from a phyloseq object.
#' @param phobj A phyloseq object containing microbiome data.
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{Genus}: Taxonomic genus name.
#'   \item \code{prevalence}: Mean prevalence of ASVs within each genus (number of samples where ASV is present).
#'   \item \code{sum_prevalence}: Total prevalence summed across all ASVs in each genus.
#' }
#' @details
#' This function first obtains the prevalence data per ASV via \code{getRelAbundanceTab}. It also calculates relative abundances internally but does not use them in the summary. The main output summarizes mean and total prevalence by genus using \code{plyr::ddply}.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   df <- getMeanRelAbundancesByGenus(phyloseq_obj)
#'   print(df)
#' }
#' }
#' @seealso
#'  \code{\link[phyloseq]{transform_sample_counts}}
#' @rdname getMeanRelAbundancesByGenus
#' @export
#' @importFrom phyloseq transform_sample_counts
#' @importFrom plyr ddply
getMeanRelAbundancesByGenus <- function(phobj){
  #Mean prevalence of ASVs by Genus (absolute prevalence, not relative)
  pre_prevalence <- getRelAbundanceTab(phobj)

  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})

  df_prevalence <- plyr::ddply(pre_prevalence, "Genus", function(df1){cbind(mean(df1$Prevalence),
                                                                      sum(df1$Prevalence))})

  names(df_prevalence) <- c("Genus", "prevalence", "Sum of prevalences")
  return(df_prevalence)
}
