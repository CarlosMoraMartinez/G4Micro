getCLRdata <- function(phobj, prevalence_lim = 0.05){
  library(CoDaSeq)
  #1) Imputate zeros
  # cmultRepl expects samples in rows and taxa in columns
  tax_matrix <- getZCompositionImputedData(phobj, prevalence_lim)

  #2) CLR transform
  tax_matrix_clr <- codaSeq.clr(tax_matrix, samples.by.row=F, IQLR=TRUE)

  return(tax_matrix_clr)
}
