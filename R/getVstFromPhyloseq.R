#' @title Variance-Stabilizing Transformation from a Phyloseq Object
#' @description Applies DESeq2's variance-stabilizing transformation (VST) to count data
#'              extracted from a \code{\link[phyloseq]{phyloseq}} object.
#'
#' @param phobj A \code{\link[phyloseq]{phyloseq}} object containing microbiome count data.
#' @param rows_as_samples Logical. If \code{TRUE} (default), the output table will have samples as rows
#'                        and features (taxa) as columns. If \code{FALSE}, the orientation is preserved
#'                        as in DESeq2 (features as rows, samples as columns).
#' @param outdir A character string specifying the output directory for the written VST-transformed counts.
#'               Default: \code{'./'}. If set to \code{""}, no file is written.
#' @param name A character string appended to the output file name. Default: \code{''}.
#'
#' @return A data frame containing variance-stabilized counts.
#'         The first column is \code{sampleID}, followed by transformed abundance values.
#'
#' @details The function converts a \code{phyloseq} object to a \code{DESeq2} dataset,
#'          estimates size factors and dispersions, and applies the variance-stabilizing
#'          transformation (VST) with \code{blind = TRUE}.
#'          The resulting VST counts are returned as a data frame and optionally written
#'          to a tab-separated file (\code{VST_counts_<name>.tsv}) in \code{outdir}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(phyloseq)
#'   data(GlobalPatterns)
#'   # Apply VST
#'   vst_res <- getVstFromPhyloseq(GlobalPatterns, rows_as_samples = TRUE, name = "GP")
#'   head(vst_res)
#' }
#' }
#'
#' @seealso
#'  \code{\link[phyloseq]{phyloseq_to_deseq2}},
#'  \code{\link[DESeq2]{varianceStabilizingTransformation}},
#'  \code{\link[DESeq2]{estimateSizeFactors}},
#'  \code{\link[DESeq2]{estimateDispersions}}
#' @rdname getVstFromPhyloseq
#' @export
#' @importFrom phyloseq phyloseq_to_deseq2
#' @importFrom DESeq2 estimateSizeFactors estimateDispersions
getVstFromPhyloseq <- function(phobj, rows_as_samples=TRUE, outdir="./", name=""){
  dds <- phyloseq::phyloseq_to_deseq2(phobj, ~ 1)  # no design needed, just for VST
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)
  vst_counts <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

  if(rows_as_samples){
    vstdf_write <- vst_counts %>% t %>% as.data.frame %>%
      rownames_to_column("sampleID")
  }else{
    vstdf_write <- vst_counts %>% as.data.frame %>%
      rownames_to_column("sampleID")
  }
  if(outdir != ""){
    write_tsv(vstdf_write, paste0(outdir, "VST_counts_", name, ".tsv"))
  }
  return(vstdf_write)
}
