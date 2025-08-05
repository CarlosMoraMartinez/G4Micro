
makeLinePlotComparingPhobjs <- function(all_model_results, opt, models_name1="padj_taxa_res", models_name2="praw_taxa_res"){

  all_sig_tables <- names(all_model_results) %>% map(\(name){
    a <- all_model_results[[name]][[models_name1]]$modummary %>%
      dplyr::mutate(taxa_group="padj")
    b <- all_model_results[[name]][[models_name2]]$modummary %>%
      dplyr::mutate(taxa_group="praw")
    rbind(a, b) %>% dplyr::mutate(input_data=name)
  }) %>% bind_rows() %>% dplyr::arrange(desc(Accuracy_l1out)) %>%
    dplyr::select(input_data, model, taxa_group, everything())

  write_tsv(all_sig_tables, file = paste0(opt$out, "/all_model_summaries.tsv"))

  (g1 <- ggplot(all_sig_tables, aes(x=model,
                                    y=Accuracy_l1out,
                                    col=input_data,
                                    fill=input_data,
                                    group=input_data))+
      facet_grid(taxa_group~.)+
      geom_point()+
      geom_line()+
      scale_color_cosmic()+
      scale_fill_cosmic() +
      ggpubr::theme_pubr() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  )
  ggsave(paste0(opt$out, "/all_model_accuracy.pdf"), g1,
         width = 8, height = 8)
  return(g1)
}
