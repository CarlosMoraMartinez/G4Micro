
#' @title Full Mediation Analysis for IMC (BMI) Effects
#' @description
#' Performs mediation analysis testing the effect of a mediator (default: IMC/BMI)
#' on the relationship between a predictor (e.g. microbial taxa) and a binary outcome (e.g. condition).
#' The function supports visualization (Venn diagrams, boxplots), power analysis,
#' and exports summary tables of mediation results.
#'
#' @param vstdf Data frame with vst-transformed expression/abundance data; contains variable names in column 'gene'.
#' @param df_all Full data frame including all variables for mediation analysis.
#' @param opt List of options, must include at least `out` (output directory path).
#' @param getVarsFunction Function that takes a summary dataframe and returns variables to test for mediation.
#' @param mediator_name Name of the mediator variable (default: "IMC").
#' @param y_name Name of the dependent/outcome variable (default: "Condition_bin").
#' @param plim p-value cutoff for selecting significant variables (default: 0.05).
#' @param plim_plot p-value cutoff for plotting (default: 0.05).
#' @param name Base name prefix for output plots and files (default: "analysis_IMC_separateModel_vjust").
#' @param wnet Width of network plots (default: 14).
#' @param hnet Height of network plots (default: 12).
#' @param wbars Width of barplots (default: 8).
#' @param hbars Height of barplots (default: 10).
#' @param wbars2 Width of secondary barplots (default: 10).
#' @param hbars2 Height of secondary barplots (default: 12).
#' @param use_color_scale Logical, whether to use color scale in plots (default: FALSE).
#' @param fix_barplot_limits Logical, whether to fix limits in barplots (default: FALSE).
#' @param custom_colors Optional vector of custom colors for plots (default: NULL).
#' @param make_boxplots Logical, whether to create boxplots (default: TRUE).
#' @param list2merge Optional named list of data.frames with p-values to merge for analysis; if NULL, defaults will be used.
#' @param make_power_test Logical, whether to perform power analysis on mediation results (default: FALSE).
#' @param min_power Minimum power threshold for power analysis (default: 0.9).
#'
#' @return
#' A list containing:
#' - plotNet: List of plots including Venn diagrams and mediation network plots.
#' - barplots: List of boxplots and barplots (if `make_boxplots` is TRUE).
#' - results: Data frame of mediation analysis results.
#' - power: Power analysis results (if `make_power_test` is TRUE).
#'
#' @details
#' The function expects the `vstdf` variable names to be cleaned and matched to variables in `df_all`.
#' It merges p-values from contrasts (default or user-provided `list2merge`),
#' applies the user-defined variable selection function (`getVarsFunction`),
#' runs simple mediation models using lavaan,
#' generates plots, saves results and optionally runs power analyses.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   res <- makeFullMediationAnalysisIMC(vstdf, df_all, opt, getVarsFunction)
#' }
#' }
#' @seealso
#' \code{\link[assertthat]{assert_that}},
#' \code{\link[dplyr]{filter}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{summarise_all}},
#' \code{\link[tidyr]{unite}},
#' \code{\link[ggplot2]{ggplot}}
#'
#' @rdname makeFullMediationAnalysisIMC
#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr filter select mutate summarise_all pull
#' @importFrom tidyr unite gather spread
makeFullMediationAnalysisIMC <- function(vstdf, df_all, opt, getVarsFunction, mediator_name="IMC", y_name="Condition_bin",
                                         plim=0.05, plim_plot=0.05, name="analysis_IMC_separateModel_vjust",
                                         wnet=14, hnet=12, wbars=8, hbars=10, wbars2=10, hbars2=12, use_color_scale=FALSE,
                                         fix_barplot_limits=FALSE, custom_colors=NULL, make_boxplots=TRUE, list2merge=NULL,
                                         make_power_test = FALSE, min_power=0.9){
  print(paste0("OUTPUT: ", opt$out))
  summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>%
                             gsub("\\[|\\]", "", .) %>%
                             gsub("\\._", "_", .) %>%
                             gsub("\\/|\\(|\\)", ".", .) %>%
                             gsub("\\.$", "", .)
  )

  assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

  if(is.null(list2merge)){
    list2merge <- list(
      depr_only_padj = dea2contrasts$firstContrast$resdf,
      imc_only_padj = dea2contrasts$contrastlist2$BMI_alone$resdf,
      depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrBMI$resdf,
      #depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
      imc_adjdepr_padj = dea2contrasts$contrastlist2$BMI_corrCond$resdf
    )
  }

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
  vars2venn <- list(
    "D vs C" = summary_df %>% dplyr::filter(depr_only_padj < plim & !is.na(depr_only_padj)) %>% pull(variable),
    "D vs C adj. BMI" = summary_df %>% dplyr::filter(depr_adjimc_padj < plim & ! is.na(depr_adjimc_padj)) %>% pull(variable),
    "BMI" = summary_df %>% dplyr::filter(imc_only_padj < plim & ! is.na(imc_only_padj)) %>% pull(variable),
    "BMI adj. Depr" = summary_df %>% dplyr::filter(imc_adjdepr_padj < plim & !is.na(imc_adjdepr_padj)) %>% pull(variable)
  )

  gv <- makeVenn(vars2venn, "VennDiagram_Cond_CondAdj_BMI_BMIAdj", opt)
  gv <- makeVenn(vars2venn[1:2], "VennDiagram_Cond_CondAdj", opt)
  gv <- makeVenn(vars2venn[c(1,3,4)], "VennDiagram_Cond_BMI_BMIAdj", opt)

  medresults <- lapply(vars2test, \(x, df_all){
    df <- df_all %>% dplyr::select(all_of(c(x, y_name, mediator_name)))
    with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
    df <- df[!with_nas, ]
    param_names <- c("b", "cp", "a", "a*b", "cp+a*b")
    res <- tryCatch({makeMediationSimpleLavaan(df, x, y_name, mediator_name)}, error=\(x){
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

  ## Make power analysis
  if(make_power_test){
    cat("--> Entering Power Analysis\n")
    medPowertest <- map2(medresults$Xvar, medresults$fullmodel, \(x, medres_x){
      #print(class(medres_x))
      df <- df_all %>% dplyr::select(all_of(c(x, y_name, mediator_name)))
      with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
      df <- df[!with_nas, ]
      #medres_x <- medresults_merged %>% filter(Xvar == x)
      #param_ests <- c(a = medres_x$a_Estimate, b=medres_x$b_Estimate, cp=medres_x$cp_Estimate)
      cat("--> Calculating power for ", x, "\n")
      res <- tryCatch({makeMediationSimplePowerCurveLavaan(df, x, y_name,
                                                           mediator_name,
                                                           med_res =medres_x,
                                                           nrep=100,
                                                           min_n=100, max_n=1000, inc_n=50, error=0.1) %>%
          dplyr::mutate(taxa = x) %>%
          dplyr::select(taxa, everything())

      }, error=\(x){
        warning(x)
        xx <- data.frame(taxa = character(0),
                         param = character(0),
                         ssize=numeric(0),
                         mean_est = character(0),
                         mean_se=numeric(0),
                         mean_zscore=character(0),
                         n = numeric(0),
                         power=numeric(0))
        return(list(pwres = NULL, pwsum=NULL, pwcurve = xx))
      })
      return(res)
    }) %>% bind_rows()
    write_tsv(medPowertest, file = paste0(opt$out, "mediation_analysis_IMC_PowerAnalysis.tsv"))

  }else{
    medPowertest <- list()
  }

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
  barplot_list <- list()
  if(make_boxplots){
    barplot_list$plots1 <- makePlotBySpecies(bacnames, df_all, opt$out,
                                             paste0("IMC_separate_BySpeciesPearson", as.character(plim_plot)),
                                             quantvar="IMC_log",
                                             quantvar_name = "log(IMC)",
                                             corrmethod = "pearson", w=wbars, h=hbars, plim_plot = plim_plot)
    barplot_list$plots2 <- makePlotBySpecies(bacnames,
                                             df_all,
                                             opt$out,
                                             paste0("IMC_separate_BySpeciesSpearman", as.character(plim_plot)),
                                             quantvar="IMC_log",
                                             quantvar_name = "log(IMC)", corrmethod = "spearman", w=wbars, h=hbars, plim_plot = plim_plot)
    barplot_list$plots3 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("IMC_separate_BySpeciesKendall", as.character(plim_plot)), quantvar="IMC_log",
                                             quantvar_name = "log(IMC)", corrmethod = "kendall", w=wbars, h=hbars, plim_plot = plim_plot)

    barplot_list$plots4 <- makePlotBySpeciesEffects_BMI(bacnames, df_all, res_final,
                                                        opt$out,
                                                        paste0("IMC_separate_BySpeciesEffects_", as.character(plim_plot)),
                                                        w=wbars2, h=hbars2, plim_plot=plim_plot, fix_limits = fix_barplot_limits)
    barplot_list$plots5 <- makePlotBySpeciesEffects_BMI2(bacnames, df_all, res_final,
                                                         opt$out,
                                                         paste0("IMC_separate_BySpeciesEffects2_", as.character(plim_plot)),
                                                         w=wbars2, h=hbars2, plim_plot=plim_plot, fix_limits = fix_barplot_limits)
  }
  return(list(plotNet=plotres, barplots = barplot_list, res_final=res_final, power_analysis=medPowertest))
}

