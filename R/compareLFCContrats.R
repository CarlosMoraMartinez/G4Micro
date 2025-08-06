#' @title Compare Log Fold Changes Across Multiple Contrasts
#' @description Visualizes and compares log2 fold changes (LFC) from multiple DESeq2 contrast results,
#' filtering by significance and ordering taxa for clear comparison.
#' @param contrastlist A named list of contrast results objects, each containing a data frame with DESeq2 results (e.g. "resdf").
#' @param firstContrast A DESeq2 contrast results object used as the reference or main contrast.
#' @param contrastNamesOrdered A character vector defining the order of contrast names in the plot.
#' @param mainContrastName Character string name for the main contrast (used for labeling).
#' @param plim_select Numeric cutoff for adjusted p-value to select taxa to plot. Default is 1e-06 (in any of the contrasts).
#' @param plim_plot Numeric cutoff for adjusted p-value to determine significance coloring in the plot. Default is 0.05.
#' @param name2remove Character string pattern to remove from contrast names in the plot legend. Default is "".
#' @param resdfname Character string name of the results data frame within each contrast list element. Default is "resdf".
#' @param outdir Directory path to save the plot PDF. Default is "./".
#' @param name Filename prefix for the saved plot PDF. Default is "LFC_compare".
#' @param w Numeric width of the saved plot in inches. Default is 12.
#' @param h Numeric height of the saved plot in inches. Default is 8.
#' @param scale_mode Character specifying facet scale mode: "fixed" or "free". Default is "fixed".
#' @return A ggplot2 object representing the comparison barplot.
#' @details
#' The function combines multiple DESeq2 contrast result data frames, filters taxa by adjusted p-value,
#' and visualizes log2 fold changes with significance coloring across contrasts.
#' Taxa names are cleaned to replace underscores with spaces for better readability.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example: assuming contrastlist and firstContrast have proper DESeq2 result dfs
#'   plot <- compareLFCContrats(contrastlist, firstContrast,
#'                             contrastNamesOrdered = c("Contrast1", "Contrast2"),
#'                             mainContrastName = "MainContrast")
#'   print(plot)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{mutate}}
#' @rdname compareLFCContrats
#' @export
#' @importFrom dplyr filter mutate
compareLFCContrats <- function(contrastlist, firstContrast,
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

  tax2plot <- alldeatables %>% filter(padj<=plim_select) %>% pull(taxon) %>% unique
  tax2plot %>% length
  taxorder <- firstContrast[[resdfname]] %>% arrange(desc(log2FoldChangeShrink)) %>%
    mutate(taxon = gsub("_", " ", taxon)) %>%
    filter(taxon %in% tax2plot) %>% pull(taxon)

  alldeatables_filt <- alldeatables %>%
    dplyr::filter(taxon %in% taxorder) %>%
    dplyr::mutate(taxon=factor(taxon, levels=taxorder),
                  comparison=gsub(name2remove, "", comparison),
                  comparison=gsub("_", " ", comparison),
                  comparison=factor(comparison,
                                    levels=contrastNamesOrdered),
                  UpOrDown = ifelse(log2FoldChangeShrink<0, "Down", "Up"),
                  UpOrDownSig = ifelse(padj<=plim_plot & !is.na(padj), UpOrDown, "NS"),
                  UpOrDownSig = factor(UpOrDownSig, levels=c("Up", "NS", "Down"))
    )

  g1<-ggplot(alldeatables_filt, aes(x=taxon, y=log2FoldChangeShrink, fill=UpOrDownSig)) +
    facet_grid(comparison ~ ., scales=scale_mode) +
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c(C_CASE, C_NS, C_CTRL))+ #c("firebrick1", "lightgray","steelblue2")
    mytheme +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10,
                                      colour = "black", angle = 0, face = "bold"))+
    theme(axis.text.x = element_text(size = 6,
                                     colour = "black", angle = 45,
                                     face = "italic", hjust=1, vjust=1))


  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)
  return(g1)
}
