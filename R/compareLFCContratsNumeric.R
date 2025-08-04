compareLFCContratsNumeric <- function(contrastlist, firstContrast,
                                      contrastNamesOrdered, mainContrastName,
                                      plim_select= 0.000001, plim_plot=0.05,
                                      name2remove = "",
                                      resdfname="resdf", outdir = "./", name="LFC_compare", w=12, h=8,
                                      scale_mode="fixed"){
  alldeatables <- map(names(contrastlist),
                      \(x) contrastlist[[x]][[resdfname]] %>% mutate(comparison=x)) %>%
    bind_rows()
  alldeatables <- rbind(alldeatables, firstContrast[[resdfname]] %>%
                          mutate(comparison=mainContrastName)) %>%
    mutate(taxon = gsub("_", " ", taxon))

  tax2plot <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select) %>% pull(taxon) %>% unique

  # 1) Calculate correlations
  mat2cor <- alldeatables %>% dplyr::filter(taxon %in% tax2plot) %>%
    dplyr::select(taxon, log2FoldChange, comparison) %>%
    spread(key=comparison, value=log2FoldChange) %>%
    column_to_rownames("taxon")
  cordf <- data.frame()
  i <- gsub(" ", "_", contrastNamesOrdered[1])
  for(j in names(mat2cor[, names(mat2cor)!=i])){
    mat <- mat2cor[, c(i, j)]
    mat <- mat[!apply(mat, MAR=1, \(x)any(is.na(x))), ]
    aux <- data.frame(
      Contrast1 = i,
      Contrast2 = j,
      pearson_cor = cor(mat[, i], mat[, j], method= "pearson"),
      pearson_pval = cor.test(mat[, i], mat[, j], method= "pearson")$p.value,
      spearman_cor = cor(mat[, i], mat[, j], method= "spearman"),
      spearman_pval = cor.test(mat[, i], mat[, j], method= "spearman")$p.value,
      kendall_cor = cor(mat[, i], mat[, j], method= "kendall"),
      kendall_pval = cor.test(mat[, i], mat[, j], method= "kendall")$p.value
    )
    cordf <- rbind(cordf, aux)
  }
  cordf <- cordf %>% dplyr::mutate(across(matches("_pval$"), list(ajd=\(x)p.adjust(x, method="BH"))))
  write_tsv(cordf, file = paste0(outdir, '/', name, "_Correlations.tsv"))

  # 2) Make Venn Diagram
  tax2plot <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select) %>% pull(taxon) %>% unique
  tax2plot_up <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select & log2FoldChangeShrink >0) %>% pull(taxon) %>% unique
  tax2plot_down <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select & log2FoldChangeShrink <0) %>% pull(taxon) %>% unique
  venplots <- map(names(mat2cor[, names(mat2cor)!=i]), \(x){
    tmp <-  alldeatables %>% dplyr::filter(comparison == x & padj < plim_select & ! is.na(padj))
    vars2venn <- list(
      "More in Depr" = tax2plot_up,
      "More in Ctrl" = tax2plot_down,
      " Up" = tmp %>% dplyr::filter(log2FoldChangeShrink > 0) %>% pull(taxon),
      " Down" = tmp %>% dplyr::filter(log2FoldChangeShrink < 0) %>% pull(taxon)
    )
    names(vars2venn)[3:4] <- paste0(x, names(vars2venn)[3:4])
    gvenn <- makeVennLocal(vars2venn, name=paste0(name, "_", x, "_VennDiagram"), outdir, w=8, h=8)
    return(list(groups=vars2venn, plot = gvenn))
  })
  names(venplots)<- names(mat2cor[, names(mat2cor)!=i])
  ## 3) Mosaic plot
  library(ggmosaic)

  df2mosaic <- alldeatables %>%
    dplyr::mutate(direction = ifelse(padj > plim_select | is.na(padj), "NS",
                                     ifelse(log2FoldChangeShrink < 0, "Down", "Up"))) %>%

    dplyr::select(taxon, comparison,direction) %>%
    tidyr::spread(key = comparison, value = direction) %>%
    dplyr::filter(taxon %in% tax2plot) %>%
    dplyr::mutate(across(-taxon, \(x)ifelse(is.na(x), "NS", x)),
                  across(-taxon, \(x)factor(x, levels=c("Down", "NS", "Up" ))),
                  across(!!sym(mainContrastName), \(x)factor(x, levels=c("Down", "NS", "Up")))
    )

  gmosplots <- map(names(df2mosaic)[!names(df2mosaic) %in% c("taxon", mainContrastName)], \(contr){
    gmos <- ggplot(data = df2mosaic) +
      geom_mosaic(aes(x = product(!!sym(mainContrastName), !!sym(contr)), fill=!!sym(mainContrastName)), na.rm=T) +
      geom_mosaic_text(aes(x = product(!!sym(mainContrastName), !!sym(contr)),
                           fill=!!sym(mainContrastName), label=after_stat(.wt)),
                       na.rm=T, as.label=T, size=6) +
      scale_fill_manual(values = c(C_CASE, C_NS, C_CTRL))+
      theme_classic() +
      mytheme +
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=10,face="bold"),
            legend.position = 'none')+
      xlab(gsub("_", " ", contr))+
      ylab(gsub("_", " ", mainContrastName))
    ggsave(filename = paste0(outdir, '/', name,'_',contr, "_MosaicPlots.pdf"), width = 4, height = 4)
    return(gmos)
  })
  return(list(correlations=cordf, vennplots=venplots, mosaicplots=gmosplots))
}
