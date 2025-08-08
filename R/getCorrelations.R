#' @title Compute Pairwise Correlations Between Two Matrices
#' @description
#' Calculates Pearson and Spearman correlations, corresponding p-values,
#' normality test p-values (Shapiro-Wilk), and linear model statistics
#' for all pairs of columns between two matrices.
#' Returns a detailed results data.frame and correlation matrices.
#'
#' @param mat1 A numeric matrix or data.frame. Each column is a variable.
#' @param mat2 A numeric matrix or data.frame. Each column is a variable.
#' @return A list containing:
#'   \item{res}{Data frame with correlation statistics for each variable pair.}
#'   \item{spcors}{Spearman correlation matrix (columns of mat1 vs columns of mat2).}
#'   \item{pcors}{Pearson correlation matrix (columns of mat1 vs columns of mat2).}
#'   \item{sp_pvals}{Matrix of Spearman correlation p-values.}
#'   \item{pr_pmat}{Matrix of Pearson correlation p-values.}
#'
#' @details
#' For each pair of variables (columns from mat1 and mat2), the function:
#' - Computes Pearson and Spearman correlations and their p-values.
#' - Tests normality of each variable with Shapiro-Wilk test.
#' - Fits a linear model and extracts R-squared and slope p-value.
#' - Adjusts p-values across all tests for multiple testing using Benjamini-Hochberg.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   mat1 <- matrix(rnorm(100*3), ncol=3)
#'   colnames(mat1) <- c("VarA1", "VarA2", "VarA3")
#'   mat2 <- matrix(rnorm(100*2), ncol=2)
#'   colnames(mat2) <- c("VarB1", "VarB2")
#'   cor_results <- getCorrelations(mat1, mat2)
#'   print(cor_results$res)
#' }
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
