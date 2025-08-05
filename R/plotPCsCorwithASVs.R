#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param tab PARAM_DESCRIPTION
#' @param contrres PARAM_DESCRIPTION
#' @param num_model PARAM_DESCRIPTION, Default: 1
#' @param component PARAM_DESCRIPTION, Default: 'Comp.2'
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
#' @rdname plotPCsCorwithASVs
#' @export 
#' @importFrom dplyr mutate
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
