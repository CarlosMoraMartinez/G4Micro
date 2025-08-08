#' @title Plot Correlations Between PCs and ASV Abundances
#' @description Generates scatterplots showing the relationship between a selected PCA component and log-transformed CLR abundances of selected ASVs,
#'   faceted by ASV and sample group, with regression lines and statistics.
#'
#' @param tab Data frame containing abundance data with columns including 'OTU', 'Abundance', 'Genus', and 'Psoriasis'.
#' @param contrres Data frame containing contrasts/results with an element 'asvs_p05' that includes ASV identifiers separated by underscores.
#' @param num_model Numeric index selecting which element of \code{contrres$asvs_p05} to use. Default is 1.
#' @param component Character string naming the PCA component column in \code{tab} to correlate against. Default is "Comp.2".
#'
#' @return A ggplot2 object showing scatterplots of the selected PCA component vs. log-transformed CLR abundance of ASVs,
#'   faceted by ASV and group, with linear regression lines and polynomial equation stats.
#'
#' @details
#' The function extracts up to the first 10 ASVs from \code{contrres$asvs_p05[num_model]} and filters \code{tab} for these ASVs.
#' Abundances are CLR-transformed (log ratio to geometric mean), and the genus name is extracted to create a combined label for faceting.
#' The plot includes points, linear regression smoothing, and R² and p-value annotations.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assuming 'tab' and 'contrres' are loaded:
#'   p <- plotPCsCorwithASVs(tab, contrres)
#'   print(p)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{mutate}}, \code{\link[ggpmisc]{stat_poly_eq}}, \code{\link[ggplot2]{ggplot}}
#'
#' @rdname plotPCsCorwithASVs
#' @export
#' @importFrom dplyr mutate
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom ggpubr theme_pubr
plotPCsCorwithASVs <- function(tab, contrres, num_model=1, component="Comp.2"){
  asvs <- strsplit(contrres$asvs_p05[num_model], "_")[[1]]
  if(length(asvs)>10) asvs <- asvs[1:10]

  tab2 <- tab %>% filter(OTU %in% asvs) %>%
    group_by(OTU) %>%
    dplyr::mutate(Abundance2= log( (Abundance + 1)/exp(mean(log(Abundance+1))) ),
                  Genus2 = sapply(Genus, FUN=function(x){y<-strsplit(x, " ")[[1]]; return(y[length(y)])} ),
                  OTU2 = paste(OTU, Genus2, Psoriasis, sep=":"))

  g0 <- ggplot(tab2 , aes_string(x=component, y="Abundance2", col="Psoriasis", fill="Psoriasis"))+
    facet_wrap(OTU2 ~ ., scales = "free", ncol=4) +
    geom_point() +
    geom_smooth(method="lm") +
    stat_poly_eq(use_label(c("R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    theme_pubr() +
    #scale_color_lancet() +
    #scale_fill_lancet() +
    ylab("CLR Abundance")


  return(g0)
}
