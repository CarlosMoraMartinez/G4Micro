makeAllPCAs <- function(phobj, counts_df, genes, vars2pca, opt, name = "PCAs"){
  design <- sample_data(phobj)
  design<- data.frame(design)
  nreads <- otu_table(phobj) %>% colSums()

  design <- design %>%
    dplyr::mutate(reads_log10_current = log10(nreads[as.character(sampleID)]))

  vars2pca <- c(vars2pca,"reads_log10_current")

  pca_plots <- lapply(vars2pca,  FUN=function(vv, counts, design, genes){
    plotPCA(counts,design, genes, vv)
  },counts_df, design, genes)
  names(pca_plots) <- vars2pca

  save(pca_plots, file = paste0(opt$out, name, "PCA.RData"))
  write_tsv(as.data.frame(pca_plots[[1]]$pca$x), file = paste0(opt$out, name, "PCAnewvars.tsv"))
  write_tsv(as.data.frame(pca_plots[[1]]$pca$rotation), file = paste0(opt$out, name, "PCArotation.tsv"))

  pdf(paste0(opt$out, name, ".pdf"), width=16, height = 12)
  for(pl in vars2pca){
    tryCatch({print(pca_plots[[pl]][["plots"]])}, error = function(x){print(getNullPlot(opt, pl))})
  }
  tmp <- dev.off()

  return(pca_plots)
}
