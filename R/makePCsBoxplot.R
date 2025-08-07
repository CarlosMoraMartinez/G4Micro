#' @title Plot Boxplots of Principal Components by Condition with Significance Annotations
#' @description
#' Generates boxplots of principal component (PC) scores grouped by a categorical condition,
#' adds significance annotations for pairwise comparisons, and saves the plot and data to files.
#'
#' This function adjusts PC score directions for easier visualization when comparing two groups (e.g., "Control" vs "Depression").
#' It re-labels PCs with explained variance percentages and outputs both the plot and the underlying data.
#'
#' @param phname Character. Name/key identifying the phyloseq object or dataset within \code{all_model_results}.
#' @param all_model_results Nested list containing model summaries, PCA results, and metadata.
#' @param opt List or similar object with at least an \code{out} element specifying the output directory.
#' @param get_pcnames_from Character. Key within \code{all_model_results[[phname]]} from which to extract PC names. Default is \code{"padj_taxa_res"}.
#' @param pca_name Character. Key within \code{all_model_results[[phname]]} where PCA objects are stored. Default is \code{"padj_taxa_pcas"}.
#' @param varname Character. Name of the metadata variable used to group samples (e.g., experimental condition). Default is \code{"Condition"}.
#' @param w Numeric. Width of the saved plot in inches. Default is 4.
#' @param h Numeric. Height of the saved plot in inches. Default is 6.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{plot}{A ggplot2 object showing PC score distributions by condition with significance annotations.}
#'   \item{tab}{A tidy data frame containing PC scores, sample IDs, and condition labels used for plotting.}
#'   \item{pcfactors}{A data frame with directional factors used to flip PC scores for visualization convenience.}
#' }
#'
#' @details
#' The function orders PCs numerically, generates descriptive PC labels with variance explained (using \code{\link{getPCnamesFromAllresults}}),
#' and adjusts PC score signs so that, if comparing "Control" vs "Depression", the scores for depression are aligned for easier interpretation.
#' Pairwise t-tests between conditions are computed and annotated on the plot using \code{ggsignif::stat_signif}.
#'
#' Output files include a PDF plot, and TSV files with the plotted data and PC direction factors.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{summarise}}
#'  \code{\link[tidyr]{gather}}, \code{\link[tidyr]{spread}}
#'  \code{\link[ggsignif]{stat_signif}}
#' @rdname makePCsBoxplot
#' @export
#' @importFrom dplyr select mutate summarise
#' @importFrom tidyr gather spread
#' @importFrom ggsignif stat_signif
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
