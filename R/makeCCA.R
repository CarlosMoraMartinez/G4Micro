#' @title Perform Canonical Correlation Analysis Between Microbiome and Metadata
#' @description
#' This function performs canonical correlation analysis (CCA) between a CLR-transformed
#' microbiome feature matrix and another data matrix (e.g., clinical variables), using
#' principal components of the microbiome matrix that explain up to a given proportion of variance.
#' It produces multiple ggplot visualizations of the canonical variates and saves them to PDF,
#' along with the CCA results in an `.RData` file.
#'
#' @param datamat A numeric matrix or data frame with samples in rows and variables in columns.
#' @param tax_matrix_clr A CLR-transformed microbiome feature matrix (taxa in rows, samples in columns).
#' @param metadata A data frame containing sample metadata. Must include a column named as in `metadatavar`.
#' @param metadatavar Character string giving the metadata variable to color plots by. Default: `"Psoriasis"`.
#' @param pcvar2retain Proportion of variance to retain when selecting microbiome PCs. Default: `0.9`.
#' @param outdir Output directory where plots and results will be saved. Default: `""` (current directory).
#' @param name Base name for output files (without extension). Default: `"CCA_from_CLR"`.
#'
#' @return A list containing:
#' \item{ccares}{The result of `stats::cancor` (canonical correlation analysis).}
#' \item{g1–g6}{Six ggplot objects corresponding to different CCA visualizations.}
#'
#' @details
#' The function first performs PCA on the CLR-transformed microbiome matrix and selects the number of PCs
#' needed to explain at least `pcvar2retain` proportion of variance.
#' These PCs are then correlated with the variables in `datamat` using canonical correlation analysis (`stats::cancor`).
#' The function generates six plots:
#' \enumerate{
#'   \item Boxplot of CC1 (microbiome) scores by `metadatavar`.
#'   \item Boxplot of CC1 (metadata) scores by `metadatavar`.
#'   \item Scatterplot of CC1 vs CC2 (metadata) with variable loadings (arrows).
#'   \item Scatterplot of CC1 vs CC2 (microbiome) with metadata variable loadings.
#'   \item Scatterplot of CC1 vs CC2 (metadata) with microbiome variable loadings.
#'   \item Scatterplot of CC1 vs CC2 (microbiome) with microbiome variable loadings.
#' }
#' All plots are saved to a single PDF file, and the CCA results are saved to an `.RData` file.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   res <- makeCCA(datamat, tax_matrix_clr, metadata,
#'                  metadatavar = "Psoriasis",
#'                  pcvar2retain = 0.9,
#'                  outdir = "results",
#'                  name = "myCCA")
#' }
#' }
#' @rdname makeCCA
#' @export
#' @importFrom stats prcomp cancor
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_point geom_segment geom_text stat_ellipse unit arrow
#' @importFrom ggrepel geom_text_repel
makeCCA <- function(datamat, tax_matrix_clr, metadata,metadatavar="Psoriasis",
                    pcvar2retain=0.9, outdir = "", name="CCA_from_CLR"){
  #library(C)
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

  g1 <- ggplot(cca_df, aes(x=!!sym(metadatavar), y=!!sym("CC1_X"), col=!!sym(metadatavar))) +
    geom_boxplot()+
    mytheme
  # scale_color_lancet()
  g2 <- ggplot(cca_df, aes(x=!!sym(metadatavar), y=!!sym("CC1_Y"), col=!!sym(metadatavar))) +
    geom_boxplot()+
    mytheme
  #scale_color_lancet()
  # Y points and Y variables
  g3 <- ggplot(cca_df, aes(x=!!sym("CC1_Y"), y=!!sym("CC2_Y"),
                           col=!!sym(metadatavar), fill=!!sym(metadatavar)))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_y, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_y,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme
  #scale_color_lancet()
  # X points and Y variables
  g4 <- ggplot(cca_df, aes(x=!!sym("CC1_X"), y=!!sym("CC2_X"),
                          col=!!sym(metadatavar), fill=!!sym(metadatavar)))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_y, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_y,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme
  #scale_color_lancet()
  # Y points and X variables
  g5 <- ggplot(cca_df, ae(x=!!sym("CC1_Y"), y=!!sym("CC2_Y"),
                          col=!!sym(metadatavar), fill=!!sym(metadatavar)))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_x, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_x,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme
  #scale_color_lancet()
  # X points and X variables
  g6 <- ggplot(cca_df, aes(x=!!sym("CC1_X"), y=!!sym("CC2_X"),
                           col=!!sym(metadatavar), fill=!!sym(metadatavar)))+
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
