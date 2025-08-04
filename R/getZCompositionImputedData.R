getZCompositionImputedData <- function(phobj, prevalence_lim = 0.05){
  library(zCompositions)
  otutab <- otu_table(phobj)
  #Umbral de prevalencia
  prevalence <- apply(otutab, MAR=1, FUN=function(x)sum(x>0)/length(x))
  otutab_filt <- otutab[prevalence > prevalence_lim ,]

  # CLR transform as in Xia, Sun & Xen book, pag 350

  #1) Imputate zeros
  # cmultRepl expects samples in rows and taxa in columns
  tax_matrix <- zCompositions::cmultRepl(X = t(otutab_filt), output="p-counts") %>% t

  return(tax_matrix)
}
