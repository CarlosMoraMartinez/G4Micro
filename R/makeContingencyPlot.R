makeContingencyPlot <- function(df, var1, var2, outdir, name, w=8, h=6){
  gmos <- ggplot(data = df) +
    geom_mosaic(aes(x = product(!!sym(var1), !!sym(var2)), fill=!!sym(var1)), na.rm=T) +
    geom_mosaic_text(aes(x = product(!!sym(var1), !!sym(var2)),
                         fill=!!sym(var1), label=after_stat(.wt)),
                     na.rm=T, as.label=T, size=6) +
    theme_classic() +
    mytheme +
    theme(axis.text.x=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          axis.text.y = element_text(size = 14,
                                     colour = "black", angle = 90, face = "bold"),
          legend.position = 'none')+
    xlab(gsub("_", " ", var2))+
    ylab(gsub("_", " ", var1))
  fname <- paste0(outdir, "/", name, "_mosaic.pdf")
  ggsave(filename = fname, gmos, width = w, height = h)
  return(gmos)
}
