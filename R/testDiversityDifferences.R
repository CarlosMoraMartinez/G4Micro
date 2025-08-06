#' @title testDiversityDifferences
#' @description Test Differences in Alpha Diversity Indices Across Grouping Variables
#' @param divtab A data frame containing alpha diversity indices (or other numeric variables) as columns
#' @param vars A character vector of column names in \code{divtab} corresponding to the alpha diversity indices
#' @param groupvars A character vector of column names in \code{divtab} representing the categorical variables
#' @param outdir ath to the output directory where the results table will be saved.
#' @param name Prefix for the output filename (TSV table), Default: 'alpha_diversity'
#' @return A data frame containing the results of the statistical tests for each combination
#' of variable and grouping variable. Includes ANOVA F-statistics and p-values, t-test and
#' Wilcoxon p-values (when applicable), Shapiro-Wilk normality test, and Bartlett and Levene
#' tests for homogeneity of variance. Also includes Benjamini-Hochberg corrected p-values
#' for t-tests and Wilcoxon tests
#' @details For each continuous variable in \code{vars} and each grouping variable in \code{groupvars},
#' the function:
#'
#' 1. Tests for overall differences using one-way ANOVA.
#' 2. If the grouping variable has exactly two levels:
#'    - Performs a Student's t-test and a Wilcoxon rank-sum test.
#' 3. Applies diagnostic tests:
#'    - Shapiro-Wilk test for normality.
#'    - Bartlett’s and Levene’s tests for homogeneity of variances.
#' 4. If the grouping variable has more than two levels:
#'    - Performs all pairwise comparisons between groups using t-tests and Wilcoxon tests
#'      via the helper function \code{getTestsForAllCombinations()}.
#'
#' The output is saved as a TSV file and also returned as a data frame.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  res <- testDiversityDifferences(
#'     divtab = alpha_div_table,
#'     vars = c("Shannon", "Chao1"),
#'     groupvars = c("Sex", "Treatment"),
#'     outdir = "results/",
#'     name = "diversity_stats"
#'   )
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{aov}},
#'  \code{\link[stats]{t.test}},
#'  \code{\link[stats]{wilcox.test}},
#'  \code{\link[stats]{shapiro.test}},
#'  \code{\link[stats]{bartlett.test}},
#'  \code{\link[car]{leveneTest}},
#'  \code{\link[stats]{p.adjust}}
#' @rdname testDiversityDifferences
#' @export
#' @importFrom car leveneTest
testDiversityDifferences <- function(divtab, vars, groupvars, outdir, name="alpha_diversity"){
  res <- data.frame()
  for (v in vars){
    if(length(unique(divtab[, v])) < 2 ) next
    for(g in groupvars){
      if(length(unique(divtab[, g])) < 2 | min(table(divtab[!is.na(divtab[, v]), g]))<2) next
      meantab <- tapply(divtab[, v], divtab[, g], mean, na.rm=T)
      if(any(is.na(meantab))) next
      divtab[, g] <- as.factor(divtab[, g])
      divtab[, v] <- as.numeric(divtab[, v])
      num_groups <- length(unique(divtab[!is.na(divtab[, g]), g]))
      form <- as.formula(paste0(v, " ~ ", g))
      aovres <- aov(form, divtab)
      tres <- tryCatch(ifelse(num_groups==2, t.test(form, divtab)$p.value, NA), error=function(x)return(NA))
      wres <- tryCatch(ifelse(num_groups==2, wilcox.test(form, divtab)$p.value, NA), error=function(x)return(NA))
      swres <- shapiro.test(divtab[, v])$p.value
      bt <-  bartlett.test(form, divtab)[[3]]
      levt <-  car::leveneTest(form, divtab)[1, 3]
      aux <- data.frame(
        variable = v,
        groups = g,
        comparison = "all",
        anova_F = summary(aovres)[[1]]["F value"][1, 1],
        anova_p = summary(aovres)[[1]]["Pr(>F)"][1, 1],
        t_test = tres,
        wilcox_test = wres,
        shapiro_normality_test = swres,
        bartlett_test = bt,
        levene_test = levt
      )
      res <- rbind(res, aux)
      if(num_groups > 2){
        aux_t <- getTestsForAllCombinations(divtab[, v], divtab[, g])
        aux2 <- data.frame(
          variable = v,
          groups = g,
          comparison = aux_t$groups_compared,
          anova_F = NA,
          anova_p = NA,
          t_test = aux_t$t_pval,
          wilcox_test = aux_t$wilcox_pval,
          shapiro_normality_test = aux_t$shapiro_test,
          bartlett_test = aux_t$bartlett_test,
          levene_test = aux_t$levene_test
        )
        res <- rbind(res, aux2)
      } ## if num_groups > 2
    } # For grouping variabbles
  } # For problem variables
  res$t_corrected <- p.adjust(res$t_test, method="BH")
  res$wilcox_corrected <- p.adjust(res$wilcox_test, method="BH")
  write_tsv(res, paste0(outdir, "/", name, "_QualitVarsTests.tsv"))
  return(res)
}
