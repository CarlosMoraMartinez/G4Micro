
makeFullMediationAnalysisIMC_Aggregate <- function(opt, getVarsFunctions, mediator_name="IMC", y_name="Condition_bin",
                                                   plim=0.05, plim_plot=0.05, name="analysis_IMC_separateModel_vjust",
                                                   wnet=14, hnet=12, wbars=8, hbars=10, wbars2=10, hbars2=12, use_color_scale=FALSE,
                                                   fix_barplot_limits=FALSE, custom_colors=NULL){
  cat(opt$out)
  summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>%
                             gsub("\\[|\\]", "", .) %>%
                             gsub("\\._", "_", .) %>%
                             gsub("\\/|\\(|\\)", ".", .) %>%
                             gsub("\\.$", "", .)
  )

  assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

  list2merge <- list(
    depr_only_padj = dea2contrasts$firstContrast$resdf,
    imc_only_padj = dea2contrasts$contrastlist2$BMI_alone$resdf,
    depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
    imc_adjdepr_padj = dea2contrasts$contrastlist2$BMI_corrCond$resdf
  )

  merged_pvals <- list2merge %>%
    lapply(\(x){
      x$taxon <- x$taxon %>%
        gsub("-", ".", .) %>%
        gsub("\\[|\\]", "", .) %>%
        gsub("\\._", "_", .) %>%
        gsub("\\/|\\(|\\)", ".", .) %>%
        gsub("\\.$", "", .)
      x <- x[match(summary_df$variable, x$taxon), ]
      #x$padj[is.na(x$padj)] <- 1
      return(x$padj)
    }) %>% bind_cols()
  summary_df <- cbind(summary_df, merged_pvals)

  vars2test <- getVarsFunction(summary_df)

  medresults <- lapply(vars2test, \(x, df_all){
    df <- df_all %>% dplyr::select(all_of(c(x, y_name, mediator_name)))
    with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
    df <- df[!with_nas, ]
    param_names <- c("b", "cp", "a", "a*b", "cp+a*b")
    res <- tryCatch({makeMediationSimple(df, x, y_name, mediator_name)}, error=\(x){
      xx <- data.frame(Estimate=rep(NA, 5),
                       S.E. = rep(NA, 5),
                       `Z-score`=rep(NA, 5),
                       p.value=rep(NA, 5))
      rownames(xx) <- c(param_names)
      return(list(estimates=xx))
    })
    res2 <- res$estimates %>% rownames_to_column("param") %>%
      dplyr::select(param, Estimate, p.value) %>%
      gather(key="var", value="value", Estimate, p.value) %>%
      filter(param %in% param_names) %>%
      tidyr::unite("tmp", param, var, sep="_") %>%
      spread(tmp, value) %>%
      mutate(Xvar = x, Yvar = y_name, Mediator=mediator_name, fullmodel=list(res$estimates))
    return(res2)
  }, df_all) %>% bind_rows()
  #all(medresults$Xvar %in% summary_df$variable)
  medresults_merged <- merge(medresults %>% dplyr::select(-fullmodel), summary_df, by.x="Xvar", by.y="variable")
  write_tsv(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.tsv"))
  save(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.RData"))

  ## Transform to plot
  def_effects <- c('a', 'b', 'cp','a*b', 'cp+a*b')
  res_trans <- map(1:nrow(medresults), \(i){
    res_i <- medresults$fullmodel[[i]]
    oldnames2 <- rownames(res_i) %>%
      sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1])
    resdf_i <- res_i %>%
      rownames_to_column("param") %>%
      mutate(x_labels = ifelse(param %in% def_effects,
                               medresults$Xvar[i],
                               gsub("V\\[|\\]", "", param)),
             param = ifelse(param %in% def_effects,
                            gsub("(b|a|cp)", paste0("\\1", as.character(i)), param, perl=T), param)
      )
  }) %>% bind_rows()
  bs_only <- res_trans %>% dplyr::filter(grepl("^b[0-9]+$", param))
  g_bs <- ggplot(bs_only,aes(x=Estimate))+geom_histogram()+mytheme+ggtitle("Histogram of b (IMC->Depr)")
  ggsave(filename = paste0(opt$out, "histogram_bs_IMC_separate.pdf"), g_bs)

  res_final <- rbind(
    bs_only %>%
      mutate_if(is.character, \(x)"b") %>%
      group_by(param, x_labels) %>%
      dplyr::summarise_all(mean) %>%
      dplyr::select(all_of(names(res_trans))),
    res_trans %>% dplyr::filter(!grepl("^b[0-9]+$", param)) %>%
      dplyr::mutate(param=gsub("b[0-9]+$", "b", param))
  )
  if(AJUST_PVALS){
    res_final$p_raw <- res_final$p.value
    res_final$p.value <- p.adjust(res_final$p.value, method = "BH")
  }
  write_tsv(res_final, file=paste0(opt$out, "mediation_analysis_IMC_separate_reformat.tsv"))
  write_tsv(bs_only, file=paste0(opt$out, "mediation_analysis_IMC_separate_bs.tsv"))

  plotres <- plotIMCMediationSimple(res_final, vars2test, opt$out, name, plim_plot = plim_plot,
                                    use_color_scale = use_color_scale, w=wnet, h=hnet, custom_colors=custom_colors)

  ## Make also boxplot
  bacnames <- plotres$bacorder$x_labels
  plots1 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("IMC_separate_BySpeciesPearson", as.character(plim_plot)), quantvar="IMC_log",
                              quantvar_name = "log(IMC)", corrmethod = "pearson", w=wbars, h=hbars, plim_plot = plim_plot)
  plots2 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("IMC_separate_BySpeciesSpearman", as.character(plim_plot)), quantvar="IMC_log",
                              quantvar_name = "log(IMC)", corrmethod = "spearman", w=wbars, h=hbars, plim_plot = plim_plot)
  plots3 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("IMC_separate_BySpeciesKendall", as.character(plim_plot)), quantvar="IMC_log",
                              quantvar_name = "log(IMC)", corrmethod = "kendall", w=wbars, h=hbars, plim_plot = plim_plot)

  plots4 <- makePlotBySpeciesEffects_BMI(bacnames, df_all, res_final, opt$out, paste0("IMC_separate_BySpeciesEffects", as.character(plim_plot)),
                                         w=wbars2, h=hbars2, plim_plot=plim_plot, fix_limits = fix_barplot_limits)

  return(list(plotNet=plotres, barplots = list(plots1, plots2, plots3, plots4), res_final=res_final))
}
