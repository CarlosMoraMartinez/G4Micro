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
