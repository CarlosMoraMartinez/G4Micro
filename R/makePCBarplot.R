
makePCBarplot <- function(phname, all_model_results, pcBoxplots, daa_all, opt,
                          get_pcnames_from="padj_taxa_res",
                          pca_name="padj_taxa_pcas",
                          varname="Condition", w=8, h=14){
  outdir <- paste0(opt$out, "/", phname)
  if(!dir.exists(outdir)) dir.create(outdir)
  predictions <- all_model_results[[phname]][[get_pcnames_from]]$models$`KNN-K=5`$preds_no_l1o
  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames
  pc_order<- gsub("PC", "", PCs) %>% as.numeric %>% order
  PCs <- PCs[pc_order]
  PCs_newnames <- getPCnamesFromAllresults(phname, all_model_results, get_pcnames_from, pca_name, varname)
  PCs_newnames <- PCs_newnames[pc_order]
  daatab <- daa_all[[phname]]$resdf
  pcfactors <- pcBoxplots[[phname]]$pcfactors %>% column_to_rownames("PC") %>% dplyr::select(factor)

  df <- all_model_results$remove_tanda2[[pca_name]][[varname]]$pca$rotation %>%
    as.data.frame() %>%
    dplyr::select(all_of(PCs)) %>%
    rownames_to_column("taxon") %>%
    dplyr::mutate(across(all_of(PCs), \(x)x*pcfactors[cur_column(), 1])) ## Multiplicar por factor para que coincida con LFC

  dfmerged <- merge(df, daatab, by="taxon", all.x=T, all.y=F) %>%
    dplyr::arrange(desc(!!sym(PCs[1]))) %>%
    dplyr::mutate(taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  taxon = factor(taxon, levels = taxon))

  newlevnames <- c( PCs_newnames, "LFC"="log2FoldChangeShrink","-10log(adj. p)"="padj")
  dflong <- dfmerged %>%
    dplyr::select(all_of(c("taxon", "padj", PCs, "log2FoldChangeShrink"))) %>%
    dplyr::mutate(padj = -10*log10(padj)) %>%
    tidyr::gather(key="variable", "value", -taxon) %>%
    dplyr::mutate(color = ifelse(variable == "padj",C_NS, ifelse(value < 0, C_CTRL, C_CASE)),
                  variable = fct_recode(variable, !!!newlevnames),
                  variable = factor(variable, levels = names(newlevnames))
    )

  gbars <- ggplot(dflong, aes(y=value, x=taxon, fill=color)) +
    facet_wrap(~ variable, nrow=1, scales = "free_x")+
    geom_col()+
    coord_flip() +
    theme_classic() +
    scale_fill_manual(values = c(C_CTRL, C_CTRL_LINK2, C_CASE))+
    theme(axis.text.y = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "italic"))+
    theme(axis.text.x = element_text(size = 10,
                                     colour = "black", angle = 0,
                                     face = "plain"))+
    theme(strip.text.x = element_text(size = 14,
                                      colour = "black", angle = 0, face = "plain")) +
    thin_barplot_lines +
    theme(legend.position="none")
  ggsave(filename = paste0(outdir, "/", get_pcnames_from, "_barplots_PCs_and_LFC.pdf"), gbars, width = w, height = h)
  write_tsv(dflong,  paste0(outdir, "/", get_pcnames_from, "_barplots_PCs_and_LFC_data.tsv"))
  return(gbars)
}
