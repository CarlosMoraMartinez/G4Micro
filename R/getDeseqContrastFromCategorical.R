#' @title Perform Differential Expression for a Categorical Contrast
#' @description Computes DESeq2 differential expression analysis for a categorical variable contrast and generates multiple shrinkage results.
#' @param dds A DESeqDataSet object.
#' @param lev_combin A character vector with 3 elements: the variable name, the reference level, and the test level.
#' @param opt A list of options including output directory and log2 fold-change threshold.
#' @param name A string to prefix output file names.
#' @return A list containing DESeq2 result objects and their corresponding data frames.
#' @details This function performs multiple shrinkage approaches using \code{lfcShrink} ("normal", "apeglm", and "ashr") and writes the results to disk.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example usage:
#'   results <- getDeseqContrastFromCategorical(dds, c("Condition", "Control", "Treated"), opt, "test1")
#' }
#' }
#' @seealso
#'   \code{\link[DESeq2]{results}},
#'   \code{\link[DESeq2]{lfcShrink}},
#'   \code{\link{defWriteDEAResults}}
#' @rdname getDeseqContrastFromCategorical
#' @export
#' @importFrom DESeq2 results lfcShrink
#' @importFrom stringr str_replace_all
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
