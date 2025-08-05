
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param tab PARAM_DESCRIPTION
#' @param contrres PARAM_DESCRIPTION
#' @param num_model PARAM_DESCRIPTION, Default: 1
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}
#' @rdname plotPCsCorwithASVs_PSO
#' @export 
#' @importFrom dplyr mutate
plotPCsCorwithASVs_PSO <- function(tab, contrres, num_model=1){
  asvs <- strsplit(contrres$asvs_p05[num_model], "_")[[1]]
  #if(length(asvs)>10) asvs <- asvs[1:10]

  tab2 <- tab %>% filter(OTU %in% asvs) %>%
    group_by(OTU) %>%
    dplyr::mutate(Abundance2= log( (Abundance + 1)/exp(mean(log(Abundance+1))) ),
                  Genus2 = sapply(Genus, FUN=function(x){y<-strsplit(x, " ")[[1]]; return(y[length(y)])} ),
                  OTU2 = paste(OTU, Genus2, Psoriasis, sep=":"))

  g0 <- ggplot(tab2 , aes_string(x="Psoriasis", y="Abundance2", col="Psoriasis"))+
    facet_wrap(OTU ~ ., scales = "free", ncol=4) +
    geom_boxplot() +
    geom_jitter(alpha=0.6) +
    #geom_signif(test="wilcox.test", na.rm=T, comparisons = c("yes", "no"),
    # step_increase=0.03,
    # tip_length = 0.01,
    # vjust=0.4,
    # color = "black"
    #) +
    theme_pubr() +
    #scale_color_lancet() +
    #scale_fill_lancet() +
    ylab("CLR Abundance")


  return(g0)
}
