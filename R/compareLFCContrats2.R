#' @title Compare Log Fold Changes Across Multiple Contrasts
#' @description Combines and visualizes shrunken log2 fold changes from multiple contrasts,
#' filtering taxa by adjusted p-value and ordering them based on a main contrast.
#' @param contrastlist A named list where each element contains a results data frame (e.g., "resdf") from a DESeq2 contrast, with a column \code{taxon}, \code{log2FoldChangeShrink}, and \code{padj}.
#' @param firstContrast A single contrast results object (list containing a \code{resdf}) used as the reference for ordering taxa.
#' @param contrastNamesOrdered Character vector giving the order of contrasts in the plot facets.
#' @param mainContrastName Name to assign to the \code{firstContrast} when plotting.
#' @param plim_select Adjusted p-value threshold for selecting taxa to display (applied to any contrast). Default: \code{1e-6}.
#' @param plim_plot Adjusted p-value threshold for significance coloring in the plot. Default: \code{0.05}.
#' @param name2remove Pattern to remove from contrast names before plotting. Default: \code{""}.
#' @param resdfname Name of the data frame within each list element that contains the results. Default: \code{"resdf"}.
#' @param outdir Directory where the PDF will be saved. Default: \code{"./"}.
#' @param name Filename prefix for the saved plot. Default: \code{"LFC_compare"}.
#' @param w Width of the saved PDF in inches. Default: \code{12}.
#' @param h Height of the saved PDF in inches. Default: \code{8}.
#' @param scale_mode Facet scale mode for ggplot: \code{"fixed"} or \code{"free"}. Default: \code{"fixed"}.
#' @return A \code{ggplot} object with the barplot comparison of log2 fold changes.
#' @details
#' The function merges results from all contrasts, selects taxa significant in at least one contrast, orders taxa by the main contrast, and creates a horizontal barplot with facets for each contrast. Significance is indicated by color.
#' Underscores in taxon and contrast names are replaced with spaces for readability.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   g <- compareLFCContrats2(contrastlist, firstContrast,
#'                            contrastNamesOrdered = c("Contrast1", "Contrast2"),
#'                            mainContrastName = "Main")
#'   print(g)
#' }
#' }
#' @seealso \code{\link[dplyr]{filter}}, \code{\link[dplyr]{mutate}}
#' @rdname compareLFCContrats2
#' @export
#' @importFrom dplyr filter mutate
compareLFCContrats2 <- function(contrastlist, firstContrast,
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
    facet_grid(. ~ comparison, scales=scale_mode) +
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c(C_CASE, C_NS,C_CTRL))+
    mytheme +
    theme_classic()+
    theme(strip.text.y = element_text(size = 8,
                                      colour = "black", angle = 0, face = "bold"))+
    theme(axis.text.y = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "italic", hjust=1, vjust=0.3))+
    coord_flip() +
    thin_barplot_lines

  ggsave(filename = paste0(outdir, "/", name, "2.pdf"), g1, width = w, height = h)
  return(g1)
}
