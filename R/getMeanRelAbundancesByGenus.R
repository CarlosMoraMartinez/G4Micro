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
