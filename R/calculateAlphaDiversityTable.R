#' @title Calculate Alpha-Diversity Indices for Each Sample
#' @description Computes several common alpha-diversity indices for each sample in a phyloseq object and exports the results to a tab-separated file.
#'
#' @param phseq_obj A \code{\link[phyloseq]{phyloseq}} object containing the microbiome data.
#' @param outdir Character string. Path to the output directory where the results table will be saved.
#' @param indices Character vector of diversity indices to calculate. Supported measures include:
#' \code{"Observed"}, \code{"Chao1"}, \code{"ACE"}, \code{"Shannon"}, \code{"Simpson"}, \code{"InvSimpson"}, and \code{"Fisher"}.
#' Default: \code{c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher")}.
#' @param name Character string. Base name for the output file (without extension). Default: \code{"diversity"}.
#'
#' @return A \code{data.frame} containing:
#' \itemize{
#'   \item All sample metadata from the phyloseq object
#'   \item One column for each calculated diversity index
#' }
#'
#' @details
#' This function uses \code{\link[phyloseq]{estimate_richness}} to calculate the specified alpha-diversity indices
#' for each sample in the phyloseq object.
#' The results are merged with the sample metadata and written to a tab-separated values (TSV) file in \code{outdir}.
#'
#' The \code{sampleID} column in the output table corresponds to the sample names in the phyloseq object.
#' Any leading \code{"X"} in sample IDs is removed to match common naming conventions.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(phyloseq)
#'   data(GlobalPatterns)
#'   result <- calculateAlphaDiversityTable(GlobalPatterns, outdir = tempdir())
#'   head(result)
#' }
#' }
#'
#' @seealso
#'  \code{\link[phyloseq]{estimate_richness}},
#'  \code{\link[phyloseq]{sample_data}},
#'  \code{\link[readr]{write_tsv}}
#'
#' @rdname calculateAlphaDiversityTable
#' @export
#' @importFrom phyloseq estimate_richness sample_data
#' @importFrom readr write_tsv
calculateAlphaDiversityTable <- function(phseq_obj, outdir,
                                         indices = c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher"),
                                         name="diversity"){
  #library(ggpmisc) #Add regression formula
  div <- estimate_richness(phseq_obj,
                           measures = indices)
  div$sampleID <- gsub("^X", "", rownames(div))
  div2 <- merge(data.frame(sample_data(phseq_obj)), div, by="sampleID")

  write_tsv(div2, file=paste0(outdir, "/", name, ".tsv") )
  return(div2)
}
