#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param variable PARAM_DESCRIPTION, Default: 'Condition'
#' @param outname PARAM_DESCRIPTION, Default: 'prevalence_by_genus.tsv'
#' @param oldlevs PARAM_DESCRIPTION, Default: c("Control", "Depression")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{summarise}}, \code{\link[dplyr]{mutate}}
#'  \code{\link[tidyr]{gather}}, \code{\link[tidyr]{spread}}
#' @rdname getRelAbundancesByGenusAndVariable
#' @export 
#' @importFrom dplyr summarise mutate
#' @importFrom tidyr gather spread
getRelAbundancesByGenusAndVariable <- function(phobj, variable="Condition",
                                               outname="prevalence_by_genus.tsv",
                                               oldlevs=c("Control", "Depression")){
  #Mean prevalence of Phyla


  aa <- psmelt(phobj) %>% data.frame
  aa$Condition <- aa[, variable]

  levs=c("no", "yes")
  aa <- psmelt(phobj) %>% data.frame
  aa$Condition <- levs[match(aa[, variable], oldlevs)]

  asvmean <- aa %>% group_by(OTU) %>%
    dplyr::summarise(Genus = unique(Genus),
                     mean_abund = mean(Abundance),
                     prevalence = sum(Abundance>0),
                     num_samples = length(unique(sampleID)),
                     rel_prevalence = prevalence/num_samples) %>%
    group_by(Genus ) %>%
    dplyr::summarise(mean_rel_prev = mean(rel_prevalence)
    )

  aa <- aa %>%
    group_by(sampleID, Condition, Genus) %>%
    dplyr::summarise(
      Total_Abundance = sum(Abundance),
      ASV_prevalence = sum(Abundance > 0),
      num_ASVs = n(),
      rel_ASV_prevalence = ASV_prevalence/num_ASVs
    ) %>%
    ungroup() %>%
    group_by(Condition, Genus) %>%
    dplyr::summarise(prevalence = sum(Total_Abundance>0),
                     mean_abundance = mean(Total_Abundance),
                     n_samples = length(unique(sampleID)),
                     rel_prevalence = prevalence/n_samples,
                     mean_asv_prevalence = mean(rel_ASV_prevalence)
    )


  bb <- aa %>% tidyr::gather("variable", "value", prevalence,
                             mean_abundance, rel_prevalence,
                             n_samples, mean_asv_prevalence) %>%
    unite("tmp" , Condition, variable) %>%
    tidyr::spread(tmp, value) %>%
    dplyr::mutate(total_n_samples = no_n_samples + yes_n_samples,
                  total_prevalence = yes_prevalence + no_prevalence,
                  total_rel_prevalence = total_prevalence/total_n_samples,
                  mean_asv_prevalence = (yes_mean_asv_prevalence*yes_n_samples + no_mean_asv_prevalence*no_n_samples)/total_n_samples)
  if(outname != "") write_tsv(bb, file = outname)

  return(bb)
}
