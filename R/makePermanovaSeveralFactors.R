makePermanovaSeveralFactors <- function(phobj,
                                        dist_method = "bray",
                                        seed = 123,
                                        modlist = c(),
                                        outname = "permanovas_mult.RData"){
  ## From https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
  library(stringi)

  meta <- sample_data(phobj)
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))

  res <- tibble()
  for(vars in modlist){
    set.seed(seed)
    form <- paste0("braydist ~ ", paste(vars, sep= " + ", collapse = " + "))
    mod1 <- adonis2(as.formula(form), data = sampledf, na.action=na.exclude)
    aux <- tibble(formula = form, model = list(mod1), explained=1 - mod1["Residual", "R2"])
    res <-rbind(res, aux)
  }
  save(res, file=outname)
  return(res)

}
