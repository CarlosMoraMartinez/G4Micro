getSummaryTablesDeseq <- function(res, opt){
  namestab <- c(paste("p < ", as.character(opt$pval), sep="", collapse="") ,
                paste("LFC > ", as.character(log2(opt$fc)), sep="", collapse="")
  )
  restab <- res %>% as.data.frame() %>%
    dplyr::mutate(a = pvalue <= opt$pval,
                  b = ifelse(log2FoldChange <= -log2(opt$fc),"less frequent",
                             ifelse( log2FoldChange >= log2(opt$fc), "more frequent", "equal"))
    ) %>%
    dplyr::select(a:b) %>%
    set_names(namestab) %>%
    table()
  rownames(restab) <- ifelse(rownames(restab) == "TRUE",
                             paste("p <= ", as.character(opt$pval), sep="", collapse=""),
                             paste("p > ", as.character(opt$pval), sep="", collapse=""))

  #restab %>% kable(caption="Number of taxons(species) per category using raw p-values")

  restab_adj <- res %>% as.data.frame() %>%
    dplyr::mutate(a = padj <= opt$pval,
                  b = ifelse(log2FoldChange <= -log2(opt$fc),"less frequent",
                             ifelse( log2FoldChange >= log2(opt$fc), "more frequent", "equal"))
    ) %>%
    dplyr::select(a:b) %>%
    set_names(namestab) %>%
    table()
  rownames(restab_adj) <- ifelse(rownames(restab_adj) == "TRUE",
                                 paste("p <= ", as.character(opt$pval), sep="", collapse=""),
                                 paste("p > ", as.character(opt$pval), sep="", collapse=""))

  #restab_adj %>% kable(caption="Number of taxons (species) per category using adjusted p-values")
  write_tsv(as.data.frame(restab), paste0(opt$out, "/num_diff_rawpval.tsv"))
  write_tsv(as.data.frame(restab_adj), paste0(opt$out, "/num_diff_adjpval.tsv"))
  return(list("restab"=restab, "restab_adj"=restab_adj))
}
