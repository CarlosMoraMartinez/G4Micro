#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param contrastlist PARAM_DESCRIPTION
#' @param firstContrast PARAM_DESCRIPTION
#' @param contrastNamesOrdered PARAM_DESCRIPTION
#' @param mainContrastName PARAM_DESCRIPTION
#' @param plim_select PARAM_DESCRIPTION, Default: 1e-06
#' @param plim_plot PARAM_DESCRIPTION, Default: 0.05
#' @param name2remove PARAM_DESCRIPTION, Default: ''
#' @param resdfname PARAM_DESCRIPTION, Default: 'resdf'
#' @param outdir PARAM_DESCRIPTION, Default: './'
#' @param name PARAM_DESCRIPTION, Default: 'LFC_compare'
#' @param w PARAM_DESCRIPTION, Default: 12
#' @param h PARAM_DESCRIPTION, Default: 8
#' @param scale_mode PARAM_DESCRIPTION, Default: 'fixed'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
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
