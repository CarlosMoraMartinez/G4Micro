getNullPlot <- function(opt, name="var", error=FALSE){
  if(!error){
    emptylab <- paste("Error plotting ", name, sep="", collapse="")
  }else{
    emptylab = "Error while plotting"
  }

  gempty <- ggplot() +
    geom_text(aes(x=1, y=1),label=emptylab, size=14)+
    theme_light() +
    theme(axis.text = element_blank(), axis.title = element_blank())
  return(gempty)
}
