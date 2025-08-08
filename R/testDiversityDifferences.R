#' @title Test Differences in Alpha Diversity Indices Across Groups
#' @description
#' Performs statistical tests to compare alpha diversity indices (or other numeric variables)
#' across one or more grouping (categorical) variables.
#'
#' @param divtab A \code{data.frame} containing alpha diversity indices (or other numeric variables) as columns,
#' and grouping variables as additional columns.
#' @param vars A character vector with the names of the numeric variables in \code{divtab} to be tested
#' (e.g., \code{c("Shannon", "Chao1")}).
#' @param groupvars A character vector with the names of the categorical grouping variables in \code{divtab}.
#' @param outdir Path to the output directory where the results table will be saved.
#' @param name Prefix for the output filename (TSV table), Default: \code{"alpha_diversity"}.
#'
#' @return
#' A \code{data.frame} with the results of the statistical tests for each numeric variable–grouping variable
#' combination. Columns include:
#' \itemize{
#'   \item \code{variable} – Name of the numeric variable tested.
#'   \item \code{groups} – Name of the grouping variable.
#'   \item \code{comparison} – "all" for overall tests, or pairwise group comparisons.
#'   \item \code{anova_F}, \code{anova_p} – One-way ANOVA F-statistic and p-value.
#'   \item \code{t_test}, \code{wilcox_test} – p-values for t-test and Wilcoxon rank-sum test (two groups only).
#'   \item \code{shapiro_normality_test} – p-value from Shapiro-Wilk normality test.
#'   \item \code{bartlett_test}, \code{levene_test} – p-values for homogeneity of variance tests.
#'   \item \code{t_corrected}, \code{wilcox_corrected} – Benjamini–Hochberg adjusted p-values.
#' }
#'
#' @details
#' For each numeric variable in \code{vars} and each categorical grouping variable in \code{groupvars}:
#' \enumerate{
#'   \item Performs one-way ANOVA to test overall differences.
#'   \item If the grouping variable has exactly two levels:
#'         \itemize{
#'           \item Performs Student's t-test.
#'           \item Performs Wilcoxon rank-sum test.
#'         }
#'   \item Checks assumptions:
#'         \itemize{
#'           \item Shapiro–Wilk test for normality.
#'           \item Bartlett’s and Levene’s tests for homogeneity of variance.
#'         }
#'   \item If the grouping variable has more than two levels:
#'         \itemize{
#'           \item Performs all pairwise t-tests and Wilcoxon tests using
#'                 \code{\link{getTestsForAllCombinations}}.
#'         }
#' }
#' The results are saved as a TSV file and also returned as a \code{data.frame}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   res <- testDiversityDifferences(
#'     divtab = alpha_div_table,
#'     vars = c("Shannon", "Chao1"),
#'     groupvars = c("Sex", "Treatment"),
#'     outdir = "results/",
#'     name = "diversity_stats"
#'   )
#' }
#' }
#'
#' @seealso
#'  \code{\link[stats]{aov}},
#'  \code{\link[stats]{t.test}},
#'  \code{\link[stats]{wilcox.test}},
#'  \code{\link[stats]{shapiro.test}},
#'  \code{\link[stats]{bartlett.test}},
#'  \code{\link[car]{leveneTest}},
#'  \code{\link[stats]{p.adjust}},
#'  \code{\link[readr]{write_tsv}},
#'  \code{\link{getTestsForAllCombinations}}
#'
#' @rdname testDiversityDifferences
#' @export
#' @importFrom car leveneTest
#' @importFrom stats aov t.test wilcox.test shapiro.test bartlett.test p.adjust
#' @importFrom readr write_tsv
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
