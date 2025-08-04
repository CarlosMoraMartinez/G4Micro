getRelAbundancesByPhylumAndVariable_beforeAfter <- function(phobj, variable="Condition",
                                                            outname="prevalence_by_phylum.tsv"){
  #Mean prevalence of Phyla

  aa <- psmelt(phobj) %>% data.frame
  aa$Condition <- aa[, variable]

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
    dplyr::mutate(total_n_samples = Before_n_samples + After_n_samples,
                  total_prevalence = Before_prevalence + After_prevalence,
                  total_rel_prevalence = total_prevalence/total_n_samples,
                  mean_asv_prevalence = (Before_mean_asv_prevalence*Before_n_samples + Before_mean_asv_prevalence*Before_n_samples)/total_n_samples)
  if(outname != "") write_tsv(bb, file = outname)

  return(bb)
}
