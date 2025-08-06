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
