#' @title Compute Relative Abundances by Genus and Condition
#' @description Computes and summarizes relative abundance and prevalence statistics by genus across groups defined by a sample metadata variable (e.g., Condition). Outputs a table summarizing genus-level abundance and prevalence metrics.
#' @param phobj A phyloseq object containing microbiome data.
#' @param variable A character string specifying the sample metadata variable used to define groups (e.g., "Condition"). Default: 'Condition'
#' @param outname Output filename for saving the resulting summary table as a TSV. Set to "" to avoid saving. Default: 'prevalence_by_genus.tsv'
#' @param oldlevs A character vector of the two levels of the grouping variable (e.g., c("Control", "Depression")). These will be internally recoded to 'no' and 'yes' for comparison. Default: c("Control", "Depression")
#' @return A data frame with genus-level statistics: abundance, prevalence, and relative prevalence for each group, as well as total values. Optionally, this data frame is written to a TSV file.
#' @details This function aggregates ASV-level abundances to genus level,
#' then calculates the number of samples in which each genus is detected (prevalence),
#' its mean abundance, and the relative prevalence within and across conditions.
#' The relative prevalence is defined as the fraction of samples in which the genus is present.
#' The function assumes the input `phobj` is a valid `phyloseq` object with taxonomy, abundance, and sample metadata included.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  df <- getRelAbundancesByGenusAndVariable(phyloseq_obj, variable="Condition")
#'   print(df)
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
