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
