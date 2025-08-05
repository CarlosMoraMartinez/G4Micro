
makeSankeyPlot <- function(nodelist, linklist, name, fname, outdir, fontsize=18){

  nodelist <- nodelist %>%
    dplyr::mutate(value  = ifelse(nchar(value) > 60, gsub("\\([^)]+\\)", "", value, perl=TRUE), value)
    )
  # From https://plotly.com/r/sankey-diagram/
  print(head(nodelist))
  fig <- plot_ly(
    type = "sankey",
    arrangement="fixed",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    valueformat = ".0f",
    valuesuffix = "TWh",
    node = list(
      label = nodelist$value,
      color = nodelist$color,
      value = nodelist$nodesize,
      x = nodelist$xpos,
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    link = list(
      source = linklist$source,
      target = linklist$target,
      value =  linklist$Size,
      label =  as.character(round(linklist$LFC, 2)),
      color =  linklist$color
    )
  )

  fig <- fig %>% layout(
    title = name,
    font = list(
      size = fontsize
    ),
    xaxis = list(showgrid = F, zeroline = F),
    yaxis = list(showgrid = F, zeroline = F)
  )

  #fig
  #tmp <- paste0(outdir, "/", fname, ".png")
  #save_image(fig, tmp)
  #use_virtualenv("r-reticulate")
  htmlwidgets::saveWidget(fig, file=paste0(outdir, fname, ".html"))
  return(fig)
}
