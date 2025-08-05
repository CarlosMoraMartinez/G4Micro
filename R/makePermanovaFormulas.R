#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param formulas PARAM_DESCRIPTION
#' @param dist_method PARAM_DESCRIPTION, Default: 'bray'
#' @param seed PARAM_DESCRIPTION, Default: 123
#' @param outname PARAM_DESCRIPTION, Default: 'permanovas_mult.tsv'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[phyloseq]{distance}}
#'  \code{\link[dplyr]{arrange}}
#' @rdname makePermanovaFormulas
#' @export 
#' @importFrom phyloseq distance
#' @importFrom dplyr arrange
makePermanovaFormulas <- function(phobj, formulas, dist_method = "bray", seed = 123,
                                  outname = "permanovas_mult.tsv"){
  ## From https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
  library(stringi)

  meta <- sample_data(phobj)
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))
  names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>%
    gsub(" ", "_", .)

  modelos <- list()
  res <- data.frame()
  for(form in formulas){
    set.seed(seed)
    mod1 <- adonis2(as.formula(form), data = sampledf, na.action=na.exclude)
    res <-rbind(res, adonis2table(mod1))
    modelos[[form]] <- mod1
  }
  res <- res %>% dplyr::arrange(P)
  res$padj <- p.adjust(res$P, method="BH")
  write_tsv(res, file=outname)
  return(list(res=res, modelos=modelos))

}
