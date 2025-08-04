makeCCA <- function(datamat, tax_matrix_clr, metadata,metadatavar="Psoriasis",
                    pcvar2retain=0.9, outdir = "", name="CCA_from_CLR"){
  library(C)
  #Assumes matrices have the same number of rows (subjects) and that are ordered in the same way
  txpca <- prcomp(tax_matrix_clr %>% t)
  txsum <- txpca %>% summary()
  pc2use <- colnames(txsum$importance)[txsum$importance["Cumulative Proportion",] > pcvar2retain][1]
  txpcdata <- txpca$rotation[, 1:(which(colnames(txpca$rotation) == pc2use))]

  datamat <- cbind(datamat, ifelse(metadata[rownames(datamat), "Psoriasis"] == "yes" , 10, 0))
  colnames(datamat)[ncol(datamat)] <- "Psoriasis"
  #ccares <- CCorA(datamat, txpcdata)
  #ccares <- cca(txpcdata, datamat)
  ccares <- cancor(txpcdata, datamat)

  CC1_X <- txpcdata %*% ccares$xcoef[, 1]
  CC1_Y <- datamat %*% ccares$ycoef[, 1]

  CC2_X <- txpcdata %*% ccares$xcoef[, 2]
  CC2_Y <- datamat %*% ccares$ycoef[, 2]

  cca_df <- data.frame(CC1_X=scale(CC1_X),
                       CC1_Y=scale(CC1_Y),
                       CC2_X=scale(CC2_X),
                       CC2_Y=scale(CC2_Y)) %>%
    rownames_to_column("sampleID")
  cca_df[, metadatavar] <- metadata[match(cca_df$sampleID, metadata$sampleID), metadatavar]
  arrows_y = data.frame(x=0,
                        xend=ccares$ycoef[,1]*2,
                        y=0,
                        yend=ccares$ycoef[,2]*2) %>%
    rownames_to_column("var2name")
  arrows_x = data.frame(x=0,
                        xend=ccares$xcoef[,1]*2,
                        y=0,
                        yend=ccares$xcoef[,2]*2) %>%
    rownames_to_column("var2name")
  arrows_x <- arrows_x[1:5,]

  g1 <- ggplot(cca_df, aes_string(x=metadatavar, y="CC1_X", col=metadatavar)) +
    geom_boxplot()+
    mytheme
  # scale_color_lancet()
  g2 <- ggplot(cca_df, aes_string(x=metadatavar, y="CC1_Y", col=metadatavar)) +
    geom_boxplot()+
    mytheme
  #scale_color_lancet()
  # Y points and Y variables
  g3 <- ggplot(cca_df, aes_string(x="CC1_Y", y="CC2_Y", col=metadatavar, fill=metadatavar))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_y, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_y,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme
  #scale_color_lancet()
  # X points and Y variables
  g4 <- ggplot(cca_df, aes_string(x="CC1_X", y="CC2_X", col=metadatavar, fill=metadatavar))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_y, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_y,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme
  #scale_color_lancet()
  # Y points and X variables
  g5 <- ggplot(cca_df, aes_string(x="CC1_Y", y="CC2_Y", col=metadatavar, fill=metadatavar))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_x, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_x,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme
  #scale_color_lancet()
  # X points and X variables
  g6 <- ggplot(cca_df, aes_string(x="CC1_X", y="CC2_X", col=metadatavar, fill=metadatavar))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_x, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_x,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme
  #scale_color_lancet()

  save(ccares, file = paste0(outdir, name, ".RData"))
  pdf(paste0(outdir, "/", name, "_plots.pdf"), width = 12, height = 12)
  g1
  g2
  g3
  g4
  g5
  g6
  dev.off()
  return(list(ccares=ccares, g1=g1, g2=g2, g3=g3, g4=g4, g5=g5, g6=g6))
}
