
makeLinePlotComparingSamePhobjModels_Cov<- function(phname, condnames,
                                                    all_model_results,
                                                    name, opt, w=8, h=12){
  outdir <- paste0(opt$out, phname)
  if(!file.exists(outdir)) dir.create(outdir)
  all_sig_tables <- map(condnames, \(name){
    all_model_results[[phname]][[name]]$modummary %>%
      dplyr::mutate(taxa_group=name, input_data=name)
  }) %>% bind_rows() %>% dplyr::arrange(desc(Accuracy_l1out)) %>%
    dplyr::select(input_data, model, taxa_group, everything())

  model_levels <- all_sig_tables %>% group_by(model) %>%
    dplyr::summarise(acc=max(Accuracy_l1out)) %>%
    dplyr::arrange(desc(acc)) %>% pull(model)
  all_sig_tables <- all_sig_tables %>% dplyr::mutate(model = factor(model, levels=model_levels))
  write_tsv(all_sig_tables, file = paste0(outdir, "/all_model_summaries.tsv"))

  (g1 <- ggplot(all_sig_tables, aes(x=model,
                                    y=Accuracy_l1out,
                                    col=input_data,
                                    fill=input_data,
                                    group=input_data))+
      #facet_grid(taxa_group~.)+
      geom_point()+
      geom_line()+
      scale_color_cosmic()+
      scale_fill_cosmic() +
      ggpubr::theme_pubr() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  )
  ggsave(paste0(outdir, "/all_model_accuracy_", name, ".pdf"), g1,
         width = w, height = h)
  return(g1)

}
