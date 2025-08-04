WriteManyPlots <- function(plotlist, name, outdir, w=12, h=7, separate=F, opt=list()){
  if(! dir.exists(outdir)){dir.create(outdir)}
  if(separate){
    for(p in names(plotlist)){
      fname = paste0(outdir, "/", name, "_", p, ".pdf")
      pdf(fname, width=w, height=g)
      tryCatch({print(plotlist[[p]])}, error = function(x){print(getNullPlot(opt, p))})
      dev.off()
    }
  }else{
    fname = paste0(outdir, "/", name, "_all.pdf")
    pdf(fname, width=w, height=h)
    for(p in names(plotlist)){
      tryCatch({print(plotlist[[p]])}, error = function(x){print(getNullPlot(opt, p))})
    }
    dev.off()
  }
}
