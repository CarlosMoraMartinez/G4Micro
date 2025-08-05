#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat1 PARAM_DESCRIPTION
#' @param mat2 PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getCorrelations
#' @export 
getCorrelations <- function(mat1, mat2){
  res <- data.frame()
  sp_pmat <-matrix(nrow=ncol(mat1), ncol=ncol(mat2))
  colnames(sp_pmat) <- colnames(mat2)
  rownames(sp_pmat) <- colnames(mat1)
  pr_pmat <-sp_pmat

  for(na in colnames(mat1)){
    for(nb in colnames(mat2)){
      ct <- cor.test(mat1[, na], mat2[, nb], method="pearson")
      ct2 <- cor.test(mat1[, na], mat2[, nb], method="spearman")
      sp1 <- shapiro.test(mat1[, na])
      sp2 <- shapiro.test(mat2[, nb])
      mod <- lm(mat1[, na] ~ mat2[, nb]) %>% summary()

      aux <- data.frame(
        var1 = na,
        var2 = nb,
        pearson_cor = ct$estimate,
        pearson_pval = ct$p.value,
        spearman_cor = ct2$estimate,
        spearman_pval = ct2$p.value,
        shapiro_v1 = sp1$p.value,
        shapiro_v2 = sp2$p.value,
        rsquared = mod$r.squared,
        slope_pval = mod$coefficients[2, 4]
      )
      res <- rbind(res, aux)
      #Store results also in matrices
      sp_pmat[na, nb] <- ct2$p.value
      pr_pmat[na, nb] <- ct$p.value
    }
  }
  res$pearson_adj_pval <- p.adjust(res$pearson_pval, method="BH")
  res$spearman_adj_pval <- p.adjust(res$spearman_pval, method="BH")

  spcors <- cor(mat1, mat2, method="spearman")
  pcors <-  cor(mat1, mat2, method="pearson")
  return(list(res=res, spcors=spcors, pcors=pcors, sp_pvals = sp_pmat, pr_pmat = pr_pmat))
}
