#' @title Perform Differential Expression for a Numerical Variable Contrast
#' @description Computes DESeq2 differential expression analysis for a numerical variable contrast and generates multiple shrinkage results.
#' @param dds A DESeqDataSet object.
#' @param nvarname A character string with the name of the numerical variable (coefficient) in the design to test.
#' @param opt A list of options including output directory and log2 fold-change threshold.
#' @param name A string to prefix output file names.
#' @return A list containing DESeq2 result objects and their corresponding data frames.
#' @details This function performs multiple shrinkage approaches using \code{lfcShrink} ("normal", "apeglm", and "ashr") and writes the results to disk.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getDeseqContrastFromNumerical
#' @export
#' @importFrom DESeq2 results lfcShrink
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
