
filterBySpeciesAndFreq <- function(df, opt){
  df2 <- df %>% dplyr::filter(grepl("\\|", Pathway))%>%
    dplyr::filter(!grepl("UNINTEGRATED|UNMAPPED|UNGROUPED", Pathway, perl=T))
  keep <- rowSums(df2 > opt$mincount) > opt$minsampleswithcount
  return(df2[keep, ])
}
