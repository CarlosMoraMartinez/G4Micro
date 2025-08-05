#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dds PARAM_DESCRIPTION
#' @param nvarname PARAM_DESCRIPTION
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
#' @rdname getDeseqContrastFromNumerical
#' @export 
getDeseqContrastFromNumerical <- function(dds, nvarname, opt, name){
  res <- results(dds, name = nvarname)
  resLFC <- tryCatch(lfcShrink(dds, coef = nvarname, type="normal", lfcThreshold = log2(opt$fc)),error=\(x)data.frame() ) #apeglm gives weird results
  resLFC_ape <- tryCatch(lfcShrink(dds, coef = nvarname, type="apeglm", lfcThreshold = log2(opt$fc)),error=\(x)data.frame() ) #apeglm gives weird results
  resLFC_ashr <- tryCatch(lfcShrink(dds, coef = nvarname, type="ashr", lfcThreshold = log2(opt$fc)),error=\(x)data.frame() ) #apeglm gives weird results
  resdf <- defWriteDEAResults(res, resLFC, opt, paste0(name, "_",nvarname, "_DAAshrinkNormal.tsv"))
  resdf_ape <- defWriteDEAResults(res, resLFC_ape, opt, paste0(name, "_", nvarname, "_DAAshrinkApe.tsv"))
  resdf_ashr <- defWriteDEAResults(res, resLFC_ashr, opt, paste0(name, "_", nvarname, "_DAAshrinkAshr.tsv"))

  return(list("res"=res,
              "resLFC"=resLFC,
              "resLFC_ape"=resLFC_ape,
              "resLFC_ashr"=resLFC_ashr,
              "resdf"=resdf,
              "resdf_ape"=resdf_ape,
              "resdf_shr"=resdf_ashr,
              "nvarname"=nvarname))
}
