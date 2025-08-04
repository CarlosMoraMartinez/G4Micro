makeConstrainedOrdinationAll <- function(phobj, opt, name="CCA", outdir= "CCA", variableList = list(c("Condition")),  dist_type="bray"){
  plots <- list()
  i = 1
  for(vars in variableList){
    if(length(vars) ==1 ){
      plots[[i]] <- makeConstrainedOrdinationSingleVar(phobj, variables =vars,  dist_type=dist_type)
    }else{
      plots[[i]] <- makeConstrainedOrdination(phobj, variables =vars,  dist_type=dist_type)
    }
    i <- i+1
  }
  names(plots)<- sapply(variableList, paste, sep="-", collapse="-")
  WriteManyPlots(plots, name, outdir, opt=opt)
  return(plots)
}
