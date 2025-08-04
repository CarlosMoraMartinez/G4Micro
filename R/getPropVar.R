getPropVar <- function(pca, PC){
  xx <- summary(pca)
  percc<- round(100*xx$importance["Proportion of Variance" , PC],2) %>% as.character
  labb <- paste(PC, " (", percc, "%)", sep="", collapse="")
  return(labb)
}
