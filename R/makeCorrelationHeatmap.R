makeCorrelationHeatmap <- function(mat1, mat2, cormethod="pearson",
                                   pval=0.05, select_asvs=c(), outdir="", name="corrHeatmap",
                                   clust=T, order_rows = c(), order_cols = c()){
  # select_asvs: plot these ASVs. Overrides pval
  # pval: plot ASVs with any correlation >  than this threshold
  cormat <- cor(mat1, mat2, method=cormethod)
  cor_pval <- expand.grid(colnames(mat1), colnames(mat2)) %>%
    rowwise() %>%
    mutate(pval = cor.test(mat1[, Var1], mat2[, Var2], method=cormethod)$p.value) %>%
    ungroup() %>%
    tidyr::spread(Var2, pval) %>%
    column_to_rownames("Var1") %>%
    as.matrix()
  cor_pval <- cor_pval[rownames(cormat), colnames(cormat)]


  if(length(select_asvs) > 0){
    cormat <- cormat[select_asvs, ]
    cor_pval <- cor_pval[select_asvs, ]
  }else if(pval < 1.0){
    select_asvs <- cor_pval %>% apply(MAR=1, function(x) any(x < pval)) %>% which() %>% names
    cormat <- cormat[select_asvs, ]
    cor_pval <- cor_pval[select_asvs, ]
  }
  cor_ptext <- getSignif2(cor_pval)
  if(! clust){
    cormat <- cormat[order_rows, order_cols]
    cor_ptext <- cor_ptext[order_rows, order_cols]
    cor_pval <- cor_pval[order_rows, order_cols]
  }

  cormat <- cormat %>% t
  cor_ptext <- cor_ptext %>% t
  cor_pval <- cor_pval %>% t

  oname <- paste0(outdir, "/", name, ".pdf")
  fontsize_row = 16 - nrow(cormat) / 15
  fontsize_col = 16 - ncol(cormat) / 15

  pheatmap(cormat, display_numbers = cor_ptext, filename=oname,
           width = 20, height = 6,
           fontsize_row = fontsize_row,
           fontsize_number = 16,
           fontsize_col = fontsize_col,
           cluster_rows = clust, cluster_cols = clust
  )
  tmp <- tryCatch(dev.off(), error=function(x){})
  hm <- pheatmap(cormat, display_numbers = cor_ptext,
                 filename=oname,
                 width = 20, height = 6,
                 fontsize_row = fontsize_row,
                 fontsize_number = 16,
                 fontsize_col = fontsize_col,
                 cluster_rows = clust, cluster_cols = clust
  )
  return(list(hm=hm, asvs=select_asvs, cormat=cormat,
              pvals=cor_pval, cor_ptext=cor_ptext))
}
