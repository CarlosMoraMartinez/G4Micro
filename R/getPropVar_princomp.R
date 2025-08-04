getPropVar_princomp <- function(pca){
  pcts <- round(100*pca$sdev^2/sum(pca$sdev^2), 2) %>% as.character
  nnames <- gsub("Comp.", "PC", names(pca$sdev))
  labb <- paste(nnames, " (", pcts, "%)", sep="")
  names(labb)<- nnames
  return(labb)
}
