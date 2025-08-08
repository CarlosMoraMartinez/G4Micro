#' @title Phylum-Level Abundance and Prevalence Summary by Group
#' @description Computes relative abundance and prevalence of phyla across a binary variable (e.g., condition). Outputs and returns a table with detailed summary statistics.
#' @param phobj A \code{phyloseq} object containing microbiome data.
#' @param variable A character string indicating the sample metadata variable to group by. Default: 'Condition'
#' @param outname Filename to save the output summary table. Set to "" to avoid saving. Default: 'prevalence_by_phylum.tsv'
#' @param oldlevs A character vector of two values representing the original levels of the grouping variable, in order. These will be relabeled as 'no' and 'yes'. Default: c("Control", "Depression")
#' @return A data frame summarizing phylum-level statistics, including abundance, prevalence, and ASV prevalence per condition and overall.
#' @details This function:
#' \itemize{
#'   \item Converts the input variable (e.g., "Condition") to binary labels "no" and "yes"
#'   \item Calculates phylum-level abundance and prevalence per sample and condition
#'   \item Computes ASV-level prevalence per phylum
#'   \item Outputs a table with summary statistics such as total prevalence, relative prevalence, and mean ASV prevalence
#' }
#' @examples
#' \dontrun{
#' if(interactive()){
#'  library(phyloseq)
#'   # Load example dataset
#'   data(GlobalPatterns)
#'
#'   # Create a binary variable in sample_data
#'   sample_data(GlobalPatterns)$MyCondition <- ifelse(sample_data(GlobalPatterns)$SampleType == "Feces", "Control", "Depression")
#'
#'   # Run the function summarizing by MyCondition
#'   result <- getRelAbundancesByPhylumAndVariable(GlobalPatterns, variable="MyCondition", outname="")
#'
#'   # View the summary table
#'   head(result)
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{summarise}}, \code{\link[dplyr]{mutate}}
#'  \code{\link[tidyr]{spread}}
#' @rdname getRelAbundancesByPhylumAndVariable
#' @export
#' @importFrom dplyr summarise mutate
#' @importFrom tidyr spread
getRelAbundancesByPhylumAndVariable <- function(phobj, variable="Condition",
                                                outname="prevalence_by_phylum.tsv",
                                                oldlevs=c("Control", "Depression")){
  #Mean prevalence of Phyla
  levs=c("no", "yes")
  aa <- psmelt(phobj) %>% data.frame
  aa$Condition <- levs[match(aa[, variable], oldlevs)]

  asvmean <- aa %>% group_by(OTU) %>%
    dplyr::summarise(Phylum = unique(Phylum),
                     mean_abund = mean(Abundance),
                     prevalence = sum(Abundance>0),
                     num_samples = n(),
                     rel_prevalence = prevalence/n()) %>%
    group_by(Phylum) %>%
    dplyr::summarise(mean_rel_prev = 100*mean(rel_prevalence)
    )

  aa <- aa %>%
    group_by(sampleID, Condition, Phylum) %>%
    dplyr::summarise(
      Total_Abundance = sum(Abundance),
      ASV_prevalence = sum(Abundance > 0),
      num_ASVs = n(),
      rel_ASV_prevalence = ASV_prevalence/num_ASVs
    ) %>%
    ungroup() %>%
    group_by(Condition, Phylum) %>%
    dplyr::summarise(prevalence = sum(Total_Abundance>0),
                     mean_abundance = mean(Total_Abundance),
                     n_samples = length(unique(sampleID)),
                     rel_prevalence = prevalence/n_samples,
                     mean_asv_prevalence = mean(rel_ASV_prevalence)
    )


  bb <- aa %>% gather("variable", "value", prevalence,
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
