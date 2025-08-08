#' @title Compare Log Fold Changes Between Contrasts
#' @description
#' Given a list of contrasts and a main contrast, this function compares log2 fold changes
#' for significantly different taxa, computes correlation statistics, generates Venn diagrams
#' of shared/unique significant taxa, and produces mosaic plots showing directionality overlaps.
#'
#' @param contrastlist A named list of contrast result objects, each containing a data frame
#'   (named in `resdfname`) with differential abundance results.
#' @param firstContrast A single contrast result object representing the main contrast.
#' @param contrastNamesOrdered Character vector of contrast names in the desired plotting order.
#' @param mainContrastName Character string naming the main contrast in plots.
#' @param plim_select Numeric p-value threshold for selecting taxa for correlation and plots. Default: 1e-06.
#' @param plim_plot Numeric p-value threshold for visualization plots. Default: 0.05.
#' @param name2remove Character string to filter out taxa whose names contain this pattern. Default: ''.
#' @param resdfname Name of the element in each contrast object containing the results data frame. Default: 'resdf'.
#' @param outdir Output directory for saved plots and result tables. Default: './'.
#' @param name Base name for output files. Default: 'LFC_compare'.
#' @param w Numeric plot width for saved plots. Default: 12.
#' @param h Numeric plot height for saved plots. Default: 8.
#' @param scale_mode Character indicating scale behavior for plots ('fixed' or 'free'). Default: 'fixed'.
#'
#' @return A list with:
#'   \item{correlations}{A data frame of Pearson, Spearman, and Kendall correlations between contrasts.}
#'   \item{vennplots}{A named list of Venn diagram plots comparing significant taxa.}
#'   \item{mosaicplots}{A named list of mosaic plots showing overlap in directionality.}
#'
#' @details
#' The function works in three stages:
#' 1. **Correlation analysis** – Identifies significantly changing taxa in the main contrast and calculates correlations with all other contrasts.
#' 2. **Venn diagrams** – For each contrast, compares sets of significantly up- and down-regulated taxa against the main contrast.
#' 3. **Mosaic plots** – Visualizes overlap in significance and directionality between the main contrast and each other contrast.
#'
#' The function writes a TSV file with correlation statistics and saves all plots to `outdir`.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   res <- compareLFCContratsNumeric(
#'     contrastlist = contrast_results_list,
#'     firstContrast = main_contrast,
#'     contrastNamesOrdered = c("contrast1", "contrast2"),
#'     mainContrastName = "contrast1"
#'   )
#' }
#' }
#'
#' @seealso
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}},
#'  \code{\link[tidyr]{spread}}
#' @rdname compareLFCContratsNumeric
#' @export
#' @importFrom dplyr filter select mutate
#' @importFrom tidyr spread
#' @importFrom ggmosaic geom_mosaic_text geom_mosaic product
compareLFCContratsNumeric <- function(contrastlist, firstContrast,
                                      contrastNamesOrdered, mainContrastName,
                                      plim_select= 0.000001, plim_plot=0.05,
                                      name2remove = "",
                                      resdfname="resdf", outdir = "./", name="LFC_compare", w=12, h=8,
                                      scale_mode="fixed"){
  alldeatables <- map(names(contrastlist),
                      \(x) contrastlist[[x]][[resdfname]] %>% mutate(comparison=x)) %>%
    bind_rows()
  alldeatables <- rbind(alldeatables, firstContrast[[resdfname]] %>%
                          mutate(comparison=mainContrastName)) %>%
    mutate(taxon = gsub("_", " ", taxon))

  tax2plot <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select) %>% pull(taxon) %>% unique

  # 1) Calculate correlations
  mat2cor <- alldeatables %>% dplyr::filter(taxon %in% tax2plot) %>%
    dplyr::select(taxon, log2FoldChange, comparison) %>%
    spread(key=comparison, value=log2FoldChange) %>%
    column_to_rownames("taxon")
  cordf <- data.frame()
  i <- gsub(" ", "_", contrastNamesOrdered[1])
  for(j in names(mat2cor[, names(mat2cor)!=i])){
    mat <- mat2cor[, c(i, j)]
    mat <- mat[!apply(mat, MAR=1, \(x)any(is.na(x))), ]
    aux <- data.frame(
      Contrast1 = i,
      Contrast2 = j,
      pearson_cor = cor(mat[, i], mat[, j], method= "pearson"),
      pearson_pval = cor.test(mat[, i], mat[, j], method= "pearson")$p.value,
      spearman_cor = cor(mat[, i], mat[, j], method= "spearman"),
      spearman_pval = cor.test(mat[, i], mat[, j], method= "spearman")$p.value,
      kendall_cor = cor(mat[, i], mat[, j], method= "kendall"),
      kendall_pval = cor.test(mat[, i], mat[, j], method= "kendall")$p.value
    )
    cordf <- rbind(cordf, aux)
  }
  cordf <- cordf %>% dplyr::mutate(across(matches("_pval$"), list(ajd=\(x)p.adjust(x, method="BH"))))
  write_tsv(cordf, file = paste0(outdir, '/', name, "_Correlations.tsv"))

  # 2) Make Venn Diagram
  tax2plot <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select) %>% pull(taxon) %>% unique
  tax2plot_up <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select & log2FoldChangeShrink >0) %>% pull(taxon) %>% unique
  tax2plot_down <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select & log2FoldChangeShrink <0) %>% pull(taxon) %>% unique
  venplots <- map(names(mat2cor[, names(mat2cor)!=i]), \(x){
    tmp <-  alldeatables %>% dplyr::filter(comparison == x & padj < plim_select & ! is.na(padj))
    vars2venn <- list(
      "More in Depr" = tax2plot_up,
      "More in Ctrl" = tax2plot_down,
      " Up" = tmp %>% dplyr::filter(log2FoldChangeShrink > 0) %>% pull(taxon),
      " Down" = tmp %>% dplyr::filter(log2FoldChangeShrink < 0) %>% pull(taxon)
    )
    names(vars2venn)[3:4] <- paste0(x, names(vars2venn)[3:4])
    gvenn <- makeVennLocal(vars2venn, name=paste0(name, "_", x, "_VennDiagram"), outdir, w=8, h=8)
    return(list(groups=vars2venn, plot = gvenn))
  })
  names(venplots)<- names(mat2cor[, names(mat2cor)!=i])

  ## 3) Mosaic plot
  df2mosaic <- alldeatables %>%
    dplyr::mutate(direction = ifelse(padj > plim_select | is.na(padj), "NS",
                                     ifelse(log2FoldChangeShrink < 0, "Down", "Up"))) %>%

    dplyr::select(taxon, comparison,direction) %>%
    tidyr::spread(key = comparison, value = direction) %>%
    dplyr::filter(taxon %in% tax2plot) %>%
    dplyr::mutate(across(-taxon, \(x)ifelse(is.na(x), "NS", x)),
                  across(-taxon, \(x)factor(x, levels=c("Down", "NS", "Up" ))),
                  across(!!sym(mainContrastName), \(x)factor(x, levels=c("Down", "NS", "Up")))
    )

  gmosplots <- map(names(df2mosaic)[!names(df2mosaic) %in% c("taxon", mainContrastName)], \(contr){
    gmos <- ggplot(data = df2mosaic) +
      geom_mosaic(aes(x = product(!!sym(mainContrastName), !!sym(contr)), fill=!!sym(mainContrastName)), na.rm=T) +
      geom_mosaic_text(aes(x = product(!!sym(mainContrastName), !!sym(contr)),
                           fill=!!sym(mainContrastName), label=after_stat(.wt)),
                       na.rm=T, as.label=T, size=6) +
      scale_fill_manual(values = c(C_CASE, C_NS, C_CTRL))+
      theme_classic() +
      mytheme +
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=10,face="bold"),
            legend.position = 'none')+
      xlab(gsub("_", " ", contr))+
      ylab(gsub("_", " ", mainContrastName))
    ggsave(filename = paste0(outdir, '/', name,'_',contr, "_MosaicPlots.pdf"), width = 4, height = 4)
    return(gmos)
  })
  return(list(correlations=cordf, vennplots=venplots, mosaicplots=gmosplots))
}
