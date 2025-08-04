makeTotalCountPlot <- function(metadata, variable, outdir, name = ""){
  g1 <- ggplot(metadata, aes_string(x=variable, y="nreads", col=variable, fill=variable)) +
    geom_boxplot(alpha=0.6, width=0.5) +
    #scale_color_lancet() +
    #scale_fill_lancet() +
    ylab("Total counts") +
    geom_jitter() + #size=3*samplesums$IMC/max(samplesums$IMC[!is.na(samplesums$IMC)])
    mytheme
  outname <- paste0(outdir,"/", name, '_', variable, ".pdf")
  ggsave(filename=outname, g1)
  return(g1)
}
