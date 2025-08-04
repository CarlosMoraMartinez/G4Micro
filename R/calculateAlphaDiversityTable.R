calculateAlphaDiversityTable <- function(phseq_obj, outdir,
                                         indices = c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher"),
                                         name="diversity"){
  library(ggpmisc) #Add regression formula
  div <- estimate_richness(phseq_obj,
                           measures = indices)
  div$sampleID <- gsub("^X", "", rownames(div))
  div2 <- merge(data.frame(sample_data(phseq_obj)), div, by="sampleID")

  write_tsv(div2, file=paste0(outdir, "/", name, ".tsv") )
  return(div2)
}
