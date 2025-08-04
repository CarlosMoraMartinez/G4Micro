defWriteDEAResults <- function(res, resLFC, opt, name){
  fname <-paste(opt$out, name, sep="/", collapse="/")
  resdf <- res %>% as.data.frame(row.names = rownames(.)) %>%
    rownames_to_column("taxon") %>%
    dplyr::mutate(log2FoldChangeShrink = resLFC$log2FoldChange,
                  lfcSE_Shrink = resLFC$lfcSE, svalue = resLFC$svalue)
  write_tsv(resdf, fname)
  return(resdf)
}
