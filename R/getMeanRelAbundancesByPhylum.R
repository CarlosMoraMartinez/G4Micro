getMeanRelAbundancesByPhylum <- function(phobj){
  #Mean prevalence (absolute, not relative) of ASVs by Phylum
  pre_prevalence <- getRelAbundanceTab(phobj)
  df_prevalence <- ddply(pre_prevalence, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                                       sum(df1$Prevalence))})

  names(df_prevalence) <- c("Phylum", "prevalence", "um of prevalences")
  return(df_prevalence)
}
