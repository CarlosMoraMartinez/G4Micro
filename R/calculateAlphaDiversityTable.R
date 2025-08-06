#' @title calculateAlphaDiversityTable
#' @description Calculates several alpha-diversity indices.
#' @param phseq_obj Phyloseq object
#' @param outdir Output directory. The results table will be saved here.
#' @param indices List of diversity indices to calculate, Default: c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher")
#' @param name PARAM_DESCRIPTION, Default: 'diversity'
#' @return A Data Frame with all the metadata and a variable for each diversity index
#' @details alculates several alpha-diversity indices for each sample in a phyloseq object,
#' using the phyloseq \code{estimate_richness()} function.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname calculateAlphaDiversityTable
#' @export
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
