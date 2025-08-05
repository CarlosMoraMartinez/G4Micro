plotAgeMediationSimple <- function(res, vars2test, outdir, outname,
                                   plim_plot = 0.05, use_color_scale=FALSE, w=14, h=6){

  C_CTRL = C_CTRL2
  edge_width_factor <- 10 ##Constants to adjust edge width. However, not used anymore
  minq <- 0.2
  rangefactor <- 10

  bacorder <- res %>% filter(grepl("^+b[0-9]+$", param)) %>% arrange(Estimate) #^cp[0-9]+$
  assertthat::assert_that(all(vars2test %in% bacorder$x_labels) & nrow(bacorder)==length(vars2test))

  vnames <- c(bacorder$x_labels, "D", "Age")
  total_age_effect <- res %>% dplyr::filter(grepl("^cp\\+", param, perl=T)) %>%
    pull(Estimate) %>%
    mean
  vertices <- data.frame(names = gsub("_", " ", vnames),
                         xpos1 = c( rep(1, length(vnames)-2), #1:(length(vnames)-2),
                                    20, #length(vnames)+5,
                                    14),#length(vnames)+5),
                         ypos1 = c(length(vnames):3,
                                   7,
                                   3 ),
                         size = c(rep(4, length(vars2test)), 10, 10),
                         totaleffect = c(bacorder$Estimate, NA, total_age_effect),
                         total_pvalue = c(bacorder$p.value, NA, NA),
                         colors = c(ifelse(bacorder$p.value < plim_plot,
                                           ifelse(bacorder$Estimate< 0, C_CTRL, C_CASE),
                                           C_NS),
                                    C_CASE, C_OTHER),
                         xoffset = ifelse(vnames == "D", 0.5, ifelse(vnames== "BMI", 1, 0))
  )


  edgenames <- expand.grid(c("a", "b"), as.character(1:length(vars2test)) ) %>%
    as.matrix %>% apply(MAR=1,paste,sep="", collapse="")
  edgenames <- c("cp", edgenames)
  assertthat::assert_that(all(edgenames %in% res$param))

  minval <- quantile(abs(res$Estimate)*edge_width_factor, minq) #to trim edge widths


  edgetab <- res %>% dplyr::filter(param %in% edgenames) %>%
    dplyr::mutate(
      from=ifelse(grepl("^b", param), gsub("_", " ", x_labels), "Age"),
      to = ifelse(grepl("^a", param), gsub("_", " ", x_labels), "D"),
      x0 = vertices$xpos1[match(from, vertices$names)],
      y0 = vertices$ypos1[match(from, vertices$names)],
      x1 = vertices$xpos1[match(to, vertices$names)],
      y1 = vertices$ypos1[match(to, vertices$names)],
      labelnames = ifelse(to=="D" & from == "Age", paste0("c'=", as.character(round(Estimate*10,2))),
                          ifelse(to=="D", paste0("b=", as.character(round(Estimate*10,2))),
                                 paste0("a=", as.character(round(Estimate*10,2))))
      ),
      color = ifelse(Estimate < 0, C_CTRL, C_CASE),
      color = ifelse(p.value > plim_plot, C_NS, color),
      linetype = ifelse(p.value > plim_plot, 2, 1),
      width = abs(Estimate)*edge_width_factor,
      width = ifelse(width <  minval, minval,
                     ifelse(width > rangefactor*minval,
                            rangefactor*minval, width))
    ) %>% dplyr::select(from, to, everything())

  edgetab$y1[edgetab$to=="D" & edgetab$from=="Age"] = edgetab$y1[edgetab$to=="D" & edgetab$from=="Age"]-0.3
  edgetab$x1[edgetab$to=="D" & edgetab$from=="Age"] = edgetab$x1[edgetab$to=="D" & edgetab$from=="Age"]-0.3

  if(use_color_scale){
    scale_link = scales::gradient_n_pal(colours=c(C_CTRL2, C_CASE_LINK),
                                        values=c(min(c(-1,min(edgetab$Estimate))),
                                                 0,
                                                 max(c(1,max(edgetab$Estimate)))))
    edgetab <- edgetab %>% dplyr::mutate(color = ifelse(p.value > plim_plot, C_NS, scale_link(Estimate)))
    vertices <- vertices %>% dplyr::mutate(colors = c(ifelse(bacorder$p.value < plim_plot,
                                                             scale_link(totaleffect),
                                                             C_NS),
                                                      C_CASE, C_OTHER))

  }

  vertices$labelnames <- mapply(vertices$names,vertices$totaleffect, FUN=function(x, val) {
    if( x == "D"){
      return(x)
    }else if(x=="Age"){
      return("Age") #paste0("Age (", as.character(round(10*val, 2)), ")")
    }else{
      return(paste0("italic('", gsub("sp ", "sp. ", x), "')"))
    }}, SIMPLIFY = T)
  (g1 <- ggplot()+
      # geom_point(data=vertices, aes(x=xpos1, y=ypos1),
      #            size=5, col=vertices$colors)+
      # geom_segment(data = edgetab, aes(x=x0, y=y0, xend=x1, yend=y1),
      #              size = edgetab$width,
      #              col = edgetab$color,
      #              linetype = edgetab$linetype,
      #              arrow=arrow(length=unit(0.01, "npc"), type="closed")) +
      geom_textsegment(data = edgetab,
                       aes(label=labelnames,x=x0, y=y0, xend=x1, yend=y1),
                       #size = edgetab$width,
                       col = edgetab$color,
                       linetype = edgetab$linetype,
                       arrow=arrow(length=unit(0.01, "npc"), type="closed")) +
      geom_label(data=vertices, aes(label=labelnames, x=xpos1, y=ypos1),
                 parse=TRUE,
                 size=vertices$size,
                 hjust=1,
                 nudge_x = vertices$xoffset,
                 col=vertices$colors, inherit.aes = T)+
      xlim(c(min(vertices$xpos1-7), max(vertices$xpos1+7)))+
      ylim(c(-1, max(vertices$ypos1)+2))+
      theme_void()
  )
  ggsave(filename = paste0(outdir, "/", name,'_p', as.character(plim_plot), ".pdf"), g1,
         width = w, height = h)
  return(list(plot=g1, bacorder=bacorder, vertices=vertices, edges=edgetab))
}
