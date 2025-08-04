makeCCAfromMatrixVegan <- function(datamat, tax_matrix_clr, pcvar2retain=0.9, outdir = "", name="CCA_from_CLR"){
  #Assumes matrices have the same number of rows (subjects) and that are ordered in the same way
  txpca <- prcomp(tax_matrix_clr %>% t)
  txsum <- txpca %>% summary()
  pc2use <- colnames(txsum$importance)[txsum$importance["Cumulative Proportion",] > pcvar2retain][1]
  txpcdata <- txpca$rotation[, 1:(which(colnames(txpca$rotation) == pc2use))]
  ccares <- CCorA(datamat, txpcdata)

  save(ccares, file = paste0(outdir, name, ".RData"))
  pdf(paste0(outdir, name, "_BiplotVegan.pdf"), width = 12, height = 12)
  biplot(ccares)
  dev.off()
  return(ccares)
}
