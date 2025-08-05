makeLinePlotComparingSamePhobjModels<- function(phname, all_model_results, opt,
                                                w=8, h=12, get_pcnames_from="padj_taxa_res", plot_extra=FALSE,
                                                plot_indiv = TRUE, plot_normal_with_smote=FALSE,
                                                filter_out = c("Ensemble2"),
                                                order_by_measure = "Accuracy_l1out", from_smote=FALSE, name="",
                                                sel_method_name="PCA DESeq"){
  outdir <- paste0(opt$out, "/", phname, "_", name, "/")
  if(!dir.exists(outdir)) dir.create(outdir)

  indivname <- "padj_taxa_res_indiv"
  lindaname <- "padj_taxa_res05linda"
  if(from_smote){
    #get_pcnames_from <- paste0(get_pcnames_from, "_SMOTE")
    indivname <- paste0(indivname, "_SMOTE")
    lindaname <- paste0(lindaname, "_SMOTE")
  }

  resph <- all_model_results[[phname]]
  pcnames <- resph[[get_pcnames_from]]$varnames
  tabs <- resph[[get_pcnames_from]]$modummary %>% dplyr::mutate(sel_method = sel_method_name, varsused = paste(pcnames, collapse="|"))
  if(lindaname %in% names(resph) & plot_extra){
    linda_pcnames <- resph[[lindaname]]$varnames
    aux <- resph[[lindaname]]$modummary %>% dplyr::mutate(sel_method = "PCA LinDA", varsused = paste(linda_pcnames, collapse="|"))
    tabs <- rbind(
      tabs,
      aux
    )
  }
  if(indivname %in% names(resph) & plot_indiv){
    tabs <- rbind(
      tabs,
      resph[[indivname]]$allmodsum
    )
  }
  if(plot_normal_with_smote){
    get_pcnames_from_smote <- paste0(get_pcnames_from, "_SMOTE")
    aux <- resph[[get_pcnames_from_smote]]$modummary %>% dplyr::mutate(sel_method = paste0(sel_method_name, " SMOTE"), varsused = paste(pcnames, collapse="|"))
    tabs <- rbind(
      tabs,
      aux
    )
  }


  write_tsv(tabs, file = paste0(outdir, phname, "_modelSummariesWithIndividualSpecies.tsv"))
  tabs2plot <- tabs %>%
    dplyr::filter(!grepl("\\+", sel_method)) %>%
    dplyr::filter(!grepl("combined 2", sel_method)) %>%
    dplyr::filter(!grepl("PC11 score", sel_method)) %>%
    dplyr::filter(!model %in% filter_out) %>%
    dplyr::mutate(model=gsub("logistic_regression", "Logistic Regr.", model),
                  sel_method = gsub("  ", " ", sel_method),
                  sel_method = gsub("DESeq", "DESeq2", sel_method)
    )
  #modorder <- tabs2plot %>% group_by(model) %>% dplyr::summarise(maxacc = max(!!sym(order_by_measure))) %>%
  #  dplyr::arrange(desc(maxacc)) %>% pull(model)

  levorder <- tabs2plot %>% group_by(sel_method) %>% dplyr::summarise(maxacc = max(!!sym(order_by_measure))) %>%
    dplyr::arrange(desc(maxacc)) %>% pull(sel_method)

  modorder <- tabs2plot %>% filter(sel_method == levorder[1]) %>%
    dplyr::arrange(desc(!!sym(order_by_measure))) %>% pull(model)

  tabs2plot <- tabs2plot %>% dplyr::mutate(model = factor(model, levels=modorder),
                                           sel_method = factor(sel_method, levels=levorder))
  (g1 <- ggplot(tabs2plot, aes(x=model,
                               y=Accuracy_l1out,
                               col=sel_method,
                               fill=sel_method,
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
      #facet_grid(sel_method~.)+
      geom_col(alpha=1, width = 0.8, position="dodge")+
      #geom_point(size=2)+
      #geom_line(alpha=1)+
      scale_color_cosmic()+
      scale_fill_cosmic() +
      ggpubr::theme_pubr() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  )
  ggsave(paste0(outdir, "/all_model_accuracy_1.pdf"), g1,
         width = w, height = w*0.5)
  tabs2plot2 <- tabs2plot %>%  dplyr::filter(!grepl("top (5|10)", sel_method, perl=T))
  g2 <- ggplot(tabs2plot2, aes(x=model,
                               y=Accuracy_l1out,
                               col=sel_method,
                               fill=sel_method,
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    #geom_col(alpha=1, width = 0.8, position="dodge")+
    geom_point(size=2)+
    geom_line(alpha=1)+
    #scale_color_cosmic()+
    #scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    ylab("Accuracy") +
    xlab("")+
    theme(legend.position="right")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_accuracy_2.pdf"), g2,
         width = 12, height = 8)
  g3 <- ggplot(tabs2plot2, aes(x=model,
                               y=Sensitivity_l1out,
                               col=sel_method,
                               fill=sel_method,
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    #geom_col(alpha=1, width = 0.8, position="dodge")+
    geom_point(size=2)+
    geom_line(alpha=1)+
    #scale_color_cosmic()+
    #scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    ylab("Sensitivity")+
    xlab("")+
    theme(legend.position="right")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_sensitivity.pdf"), g3,
         width = 12, height = 8)

  g4 <- ggplot(tabs2plot2, aes(x=model,
                               y=Specificity_l1out,
                               col=sel_method,
                               fill=sel_method,
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    #geom_col(alpha=1, width = 0.8, position="dodge")+
    geom_point(size=2)+
    geom_line(alpha=1)+
    #scale_color_cosmic()+
    #scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    ylab("Specificity")+
    xlab("") +
    theme(legend.position="right") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_Specificity.pdf"), g4,
         width = 12, height = 8)

  g5 <- ggplot(tabs2plot2, aes(x=model,
                               y=Kappa_l1out,
                               col=sel_method,
                               fill=sel_method,
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    #geom_col(alpha=1, width = 0.8, position="dodge")+
    geom_point(size=2)+
    geom_line(alpha=1)+
    #scale_color_cosmic()+
    #scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylab("Kappa")+
    xlab("") +
    scale_y_continuous(n.breaks = 6)+
    theme(legend.position="right") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_Kappa.pdf"), g5,
         width = 6, height = 3)


  g6 <- ggplot(tabs2plot2, aes(x=model,
                               y=BalancedAccuracy_l1out,
                               col=sel_method,
                               fill=sel_method,
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    #geom_col(alpha=1, width = 0.8, position="dodge")+
    geom_point(size=2)+
    geom_line(alpha=1)+
    #scale_color_cosmic()+
    #scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylab("Balanced Accuracy")+
    xlab("") +
    scale_y_continuous(n.breaks = 6)+
    theme(legend.position="right") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_BalancedAccuracy.pdf"), g6,
         width = 6, height = 3)


  g7 <- ggplot(tabs2plot2, aes(x=model,
                               y=AUC_l1out,
                               col=sel_method,
                               fill=sel_method,
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    #geom_col(alpha=1, width = 0.8, position="dodge")+
    geom_point(size=2)+
    geom_line(alpha=1)+
    #scale_color_cosmic()+
    #scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylab("AUC")+
    xlab("") +
    scale_y_continuous(n.breaks = 6)+
    theme(legend.position="right") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_AUC.pdf"), g7,
         width = 6, height = 3)

  cw <- cowplot::plot_grid(plotlist=list(g2, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname, ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined.pdf"), width = w, height = h)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g2, g5), ncol = 1)
  pdf(paste0(outdir, "/", phname, "_all_model_combined2.pdf"), width = w, height = w)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g2, g5, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname,ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined3.pdf"), width = w, height = w*1.7)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g7, g5, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname,ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined5.pdf"), width = w, height = w*1.7)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g6, g5, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname,ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined6.pdf"), width = w, height = w*1.7)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g5 + theme(legend.position = "none"),
                                         g6 + theme(legend.position = "none"),
                                         g3 + theme(legend.position = "none"),
                                         g4 + theme(legend.position = "none"),
                                         g7 + theme(legend.position = "none")), ncol = 2)
  pdf(paste0(outdir, "/", phname,ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined4.pdf"), width = w*1.7, height = w*1.5)
  print(cw)
  dev.off()

}
