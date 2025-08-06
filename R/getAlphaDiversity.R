#' @title Plot and Analyze Alpha Diversity
#' #' @description
#' This function computes and visualizes alpha diversity indices from a phyloseq object,
#' testing for significant differences across qualitative and quantitative metadata variables.
#' It generates boxplots or scatter plots for each variable and saves the results as figures.
#'
#' @param phseq_obj A \code{phyloseq} object containing microbiome data.
#' @param vars A character vector with the names of qualitative variables (factors) to compare.
#' @param qvars A character vector with names of quantitative variables to test via regression. Default: \code{c()}
#' @param opt A list of options. Must contain at least the output directory as \code{opt$out}. Default: \code{list()}
#' @param indices A character vector with the diversity indices to calculate. Default: \code{c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher")}
#' @param name A string used as a prefix for saved plots and output files. Default: \code{"AlphaDiversity"}
#' @param signif_levels Named numeric vector specifying thresholds for significance stars (e.g., "***", "**", "*", "ns"). Default: \code{c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1)}
#' @param correct_pvalues Logical; if \code{TRUE}, p-values are Bonferroni corrected across all comparisons. Default: \code{TRUE}
#' @param correct_pvalues_indices Logical; if \code{TRUE}, additional correction is applied across diversity indices. Default: \code{FALSE}
#' @param test2show A string specifying the test function to use for group comparisons (e.g., \code{"wilcox.test"}). Default: \code{"wilcox.test"}
#' @param w Width of the saved plots in inches. Default: \code{8}
#' @param h Height of the saved plots in inches. Default: \code{4}
#'
#' @return A named list of ggplot2 objects corresponding to each plotted variable.
#' Plots are also saved to disk in the specified output folder.
#'
#' @details
#' This function handles both categorical and continuous metadata variables to assess their
#' relationship with alpha diversity. Boxplots and statistical tests (e.g., Wilcoxon) are used
#' for categorical variables, while linear regressions and scatter plots are used for continuous ones.
#' The function also applies optional multiple testing corrections and automatically filters
#' variables with only one level.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plots <- getAlphaDiversity(phseq_obj, vars = c("Group", "Sex"), qvars = c("Age"),
#'                               opt = list(out = "results/"))
#'   plots$Group  # Access plot for "Group"
#'  }
#' }
#' @seealso
#'  \code{\link[base]{subset}}
#'  \code{\link[phyloseq]{prune_samples}}
#'  \code{\link[ggsignif]{stat_signif}}
#'  \code{\link[dplyr]{mutate}}
#'  \code{\link[ggpmisc]{stat_poly_eq}},
#'  \code{\link[phyloseq]{plot_richness}}
#' @rdname getAlphaDiversity
#' @export
#' @importFrom phyloseq prune_samples
#' @importFrom ggsignif stat_signif
#' @importFrom ggpmisc stat_poly_eq use_label
#' @importFrom ggpubr theme_pubclean
getAlphaDiversity <- function(phseq_obj, vars, qvars= c(),
                              opt = list(),
                              indices=c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher"),
                              name="AlphaDiversity",
                              signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                              correct_pvalues = TRUE, correct_pvalues_indices = FALSE,
                              test2show="wilcox.test", w=8, h=4){
  outdir <- paste0(opt$out, "AlphadivPlots")
  if(! dir.exists(outdir)){dir.create(outdir)}

  #Get Alpha Diversity values
  divtab <- calculateAlphaDiversityTable(phseq_obj, outdir, indices, name)
  vars <-  map_vec(divtab[, vars], \(x)length(unique(x[!is.na(x)]))) %>%
    subset(. > 1) %>% names
  # Get statistical tests (not used later)
  alphadif <- testDiversityDifferences(divtab, indices, vars, outdir, name)
  plots <- list()
  sign_df <- tibble()
  for(v in vars){
    #Si no lo convertimos a caracter, geom_signif/stat_signif fallan
    standa1 <- sample_data(phseq_obj) %>% data.frame
    standa1 <- standa1$sampleID[!is.na(standa1[, v])]
    phseq_obj_filt <- phyloseq::prune_samples(standa1, phseq_obj)
    sample_data(phseq_obj_filt)[, v] <- sample_data(phseq_obj_filt)[, v] %>% unlist %>% as.character

    divtab2 <- divtab[!is.na(divtab[, v]), ]
    divtab2[, v] <- as.character(divtab2[, v])

    comp <- combn(unique(divtab2[, v]), 2, simplify = F)
    num_comparisons <- alphadif %>%
      filter(comparison != "all" &
               groups == v &
               variable == "Observed") %>%
      nrow
    if(correct_pvalues & num_comparisons*length(indices)>1){
      signif_levels_bonferroni <- c(signif_levels[1:3]/(num_comparisons*length(indices)), signif_levels[4])
    }else{
      signif_levels_bonferroni <- signif_levels
    }
    if(correct_pvalues_indices ){
      signif_levels_bonferroni <- c(signif_levels_bonferroni[1:3]/length(indices),
                                    signif_levels_bonferroni[4])
    }

    plots[[v]] <- plot_richness(phseq_obj_filt, x = v,
                                color = v,
                                measures = indices) +
      #ggplot(divtab, aes(x=Psoriasis, y=Shannon, col=Psoriasis, fill=Psoriasis)) +
      #facet_wrap(. ~ Sex, scales = "free") +
      geom_boxplot(aes_string(fill = v), alpha = 0.7, width=0.5) +
      #scale_color_manual(values = c("#ffafcc", "#90DBF4")) +
      #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +
      #scale_color_lancet() +
      #scale_fill_lancet() +
      labs(title = v, x = '') +
      ggpubr::theme_pubclean() +
      mytheme +
      ggsignif::stat_signif(test=test2show, na.rm=T, comparisons = comp,
                            step_increase=0.06,
                            tip_length = 0.01,
                            map_signif_level=signif_levels_bonferroni,
                            vjust=0.3,
                            color = "black"
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_blank())


    #theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
  } #Plots qualitative variables

  regressions <- testDiversityWithQuantVars(divtab, indices, qvars, outdir, paste0(name, "QuantVarsRegression"))

  for(v in qvars){
    auxtext <- regressions %>% filter(predictor == v) %>%
      dplyr::mutate(text = paste0("R^2=", as.character(round(r.squared, 2)), ", p=", as.character(round(p.value, 3)) ))

    plots[[v]] <- plot_richness(phseq_obj, x = v,
                                color = v,
                                measures = c("Observed", "Chao1", "Shannon", "InvSimpson")) +
      geom_point(aes_string(fill = v), alpha = 0.7) +
      ggpmisc::stat_poly_eq(use_label(c("R2", "p")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99)+
      #stat_poly_line(method = "lm") +
      geom_smooth(method="lm", fullrange = TRUE, linetype=1) +
      #scale_color_manual(values = c("#ffafcc", "#90DBF4")) +
      #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +
      labs(title = v, x = '') +
      #geom_text(data=auxtext, aes(y=3, x = 20, labels=text)) +

      ggpubr::theme_pubclean() +
      mytheme
    theme(plot.title = element_text(hjust = 0.5)) +
      #theme(axis.text.x = element_blank())
      theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
  }
  names(plots) <- gsub("_", "-", names(plots))
  WriteManyPlots(plots, name, outdir, w=w, h=h, separate=F)
  return(plots)
}
