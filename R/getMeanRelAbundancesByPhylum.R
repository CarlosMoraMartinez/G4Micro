#' @title Calculate Mean Prevalence of ASVs by Phylum
#' @description Computes the mean and total prevalence of Amplicon Sequence Variants (ASVs) grouped by Phylum from a phyloseq object.
#' @param phobj A phyloseq object containing microbiome data.
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{Phylum}: Taxonomic phylum name.
#'   \item \code{prevalence}: Mean prevalence of ASVs within each phylum (number of samples where ASV is present).
#'   \item \code{sum_prevalence}: Total prevalence summed across all ASVs in each phylum.
#' }
#' @details
#' This function uses \code{getRelAbundanceTab} to obtain prevalence data per ASV, then summarizes the mean and sum of prevalence by phylum using \code{plyr::ddply}.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   df <- getMeanRelAbundancesByPhylum(phyloseq_obj)
#'   print(df)
#' }
#' }
#' @rdname getMeanRelAbundancesByPhylum
#' @export
#' @importFrom plyr ddply
getMeanRelAbundancesByPhylum <- function(phobj){
  #Mean prevalence (absolute, not relative) of ASVs by Phylum
  pre_prevalence <- getRelAbundanceTab(phobj)
  df_prevalence <- plyr::ddply(pre_prevalence, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                                       sum(df1$Prevalence))})

  names(df_prevalence) <- c("Phylum", "prevalence", "um of prevalences")
  return(df_prevalence)
}
