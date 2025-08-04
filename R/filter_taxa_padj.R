filter_taxa_padj <- function(resdf, plim=0.05, fc=1){
  taxa <- resdf %>%
    dplyr::filter(padj <= plim & abs(log2FoldChangeShrink) >= log2(fc) ) %>%
    pull(taxon)
  return(taxa)
}
