makePCsBoxplot <- function(phname, all_model_results, opt,
                           get_pcnames_from="padj_taxa_res",
                           pca_name="padj_taxa_pcas",
                           varname="Condition", w=4, h=6){
  outdir <- paste0(opt$out, "/", phname)
  if(!dir.exists(outdir)) dir.create(outdir)

  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames
  pc_order<- gsub("PC", "", PCs) %>% as.numeric %>% order
  PCs <- PCs[pc_order]
  PCs_newnames <- getPCnamesFromAllresults(phname, all_model_results, get_pcnames_from, pca_name, varname)
  PCs_newnames <- PCs_newnames[pc_order]
  metadata <- all_model_results[[phname]]$metadata
  dfsamples <- all_model_results$remove_tanda2[[pca_name]][[varname]]$pca$x %>%
    as.data.frame() %>%
    dplyr::select(all_of(PCs)) %>%
    rownames_to_column("sampleID") %>%
    dplyr::mutate(Condition = metadata[[varname]][match(sampleID, metadata$sampleID)]) %>%
    tidyr::gather("PC", "score", -sampleID, -Condition) %>%
    group_by(PC)

  # Multiplicar por -1 si el componente es menor en deprimidos, para plotear más fácil

  pcfactors <- dfsamples %>% group_by(PC, Condition) %>% dplyr::summarise(media = mean(score)) %>%
    tidyr::spread(Condition, media)
  if(all(as.character(unique(dfsamples$Condition)) %in% c("Control", "Depression"))){
    pcfactors <- pcfactors %>% dplyr::mutate(factor = ifelse(Depression < Control, -1, 1))
  }else{
    pcfactors <- pcfactors %>% dplyr::mutate(factor = 1)
  }

  dfsamples <- dfsamples %>%
    dplyr::mutate(score = score*pcfactors$factor[match(PC, pcfactors$PC)]) %>%
    dplyr::mutate(PC = factor(PC, levels = PCs),
                  PC = fct_recode(PC, !!!PCs_newnames),
                  Condition = gsub("Depression", "Depr.", Condition))

  comp <- combn(unique(dfsamples$Condition), 2, simplify = F)
  signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1)

  gpcbox <- ggplot(dfsamples, aes(x=Condition, y=score)) +
    facet_grid(~PC)+
    geom_violin(aes(fill=Condition)) +
    geom_boxplot(width=0.2)+
    ggsignif::stat_signif(test="t.test", na.rm=T, comparisons = comp,
                          step_increase=0.03,
                          tip_length = 0.01,
                          map_signif_level=signif_levels,
                          vjust=0.4,
                          color = "black"
    )+
    mytheme +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(filename = paste0(outdir, "/", get_pcnames_from, "_boxplot_signif_PCs.pdf"), gpcbox, width = w, height = h)
  write_tsv(dfsamples, file = paste0(outdir, "/", get_pcnames_from, "_boxplot_signif_PCs_data.tsv"))
  write_tsv(pcfactors, file = paste0(outdir, "/", get_pcnames_from, "_boxplot_signif_PCs_PCFactors.tsv"))
  return(list(plot=gpcbox, tab=dfsamples, pcfactors=pcfactors))
}
