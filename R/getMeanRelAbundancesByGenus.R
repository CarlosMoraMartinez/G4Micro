getMeanRelAbundancesByGenus <- function(phobj){
  #Mean prevalence of ASVs by Genus (absolute prevalence, not relative)
  pre_prevalence <- getRelAbundanceTab(phobj)

  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})

  df_prevalence <- ddply(pre_prevalence, "Genus", function(df1){cbind(mean(df1$Prevalence),
                                                                      sum(df1$Prevalence))})

  names(df_prevalence) <- c("Genus", "prevalence", "Sum of prevalences")
  return(df_prevalence)
}
