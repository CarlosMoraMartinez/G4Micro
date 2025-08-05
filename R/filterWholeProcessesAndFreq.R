filterWholeProcessesAndFreq <- function(df, opt){
  df2 <- df %>% dplyr::filter(!grepl("\\|", Pathway)) %>%
    dplyr::filter(!grepl("UNINTEGRATED|UNMAPPED|UNGROUPED", Pathway, perl=T))
  keep <- rowSums(df2 > opt$mincount) > opt$minsampleswithcount
  df2 <- df2[keep, ]
  return(df2)
}
