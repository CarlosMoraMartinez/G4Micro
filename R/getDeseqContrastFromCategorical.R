#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dds PARAM_DESCRIPTION
#' @param lev_combin PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getDeseqContrastFromCategorical
#' @export 
getDeseqContrastFromCategorical <- function(dds, lev_combin, opt, name){
  contrastvec <- c(lev_combin[1], lev_combin[3], lev_combin[2])
  contrast_name <- paste0(lev_combin[1],'_', lev_combin[3], '_vs_', lev_combin[2]) %>% gsub(" ", ".", .)
  res <- results(dds, contrast = contrastvec)
  resLFC <- tryCatch(lfcShrink(dds, contrast = contrastvec, type="normal", lfcThreshold = log2(opt$fc)),error=\(x)data.frame() ) #apeglm gives weird results
  resLFC_ape <- tryCatch(lfcShrink(dds, coef = contrast_name, type="apeglm", lfcThreshold = log2(opt$fc)),error=\(x)data.frame() ) #apeglm gives weird results
  resLFC_ashr <- lfcShrink(dds, contrast = contrastvec, type="ashr", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resdf <- defWriteDEAResults(res, resLFC, opt, paste0(name, "_",contrast_name, "_DAAshrinkNormal.tsv"))
  resdf_ape <- defWriteDEAResults(res, resLFC_ape, opt, paste0(name, "_", contrast_name, "_DAAshrinkApe.tsv"))
  resdf_ashr <- defWriteDEAResults(res, resLFC_ashr, opt, paste0(name, "_", contrast_name, "_DAAshrinkAshr.tsv"))

  return(list("res"=res,
              "resLFC"=resLFC,
              "resLFC_ape"=resLFC_ape,
              "resLFC_ashr"=resLFC_ashr,
              "resdf"=resdf,
              "resdf_ape"=resdf_ape,
              "resdf_shr"=resdf_ashr,
              contrast_vec = contrastvec,
              contrast_name = contrast_name))
}
