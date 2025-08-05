#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param dist_method PARAM_DESCRIPTION, Default: 'bray'
#' @param exclude_vars PARAM_DESCRIPTION, Default: c("sampleID")
#' @param outname PARAM_DESCRIPTION
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
#' @rdname makeBetaDispTests
#' @export 
#' @importFrom phyloseq distance
makeBetaDispTests<- function(phobj, dist_method = "bray",
                             exclude_vars = c("sampleID"),
                             outname){

  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))
  vars2test <- names(sampledf)[! names(sampledf) %in% exclude_vars]
  res <- data.frame()
  for(var in vars2test){
    beta <- betadisper(braydist, sampledf$Psoriasis)
    betaper <- permutest(beta)
    aux <- data.frame(variable = var,
                      F = betaper$tab$F[1],
                      p = betaper$tab$`Pr(>F)`[1])
    res <- rbind(res, aux)
  }
  write_tsv(res, file=outname)
  return(res)
}
