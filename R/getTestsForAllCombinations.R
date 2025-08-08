#' @title Perform Multiple Statistical Tests for All Pairwise Group Comparisons
#' @description
#' Given a numeric variable and a grouping factor, this function performs a series of statistical tests
#' (t-test, Wilcoxon rank-sum test, Shapiro-Wilk normality test, Bartlett's variance test, and Levene's test)
#' for all unique pairwise group combinations.
#'
#' @param var A numeric vector containing the values to be tested.
#' @param group A factor or character vector indicating group membership for each value in \code{var}.
#' @param sep Character string used to separate group names in the \code{groups_compared} column.
#' Default: \code{"-"}.
#' @param paired Logical. Should paired tests be performed? Default: \code{FALSE}.
#'
#' @return A \code{data.frame} where each row corresponds to a unique pairwise group comparison and contains:
#' \itemize{
#'   \item \code{groups_compared} — The two groups being compared, separated by \code{sep}
#'   \item \code{t_pval} — p-value from a two-sample t-test (paired or unpaired)
#'   \item \code{wilcox_pval} — p-value from a Wilcoxon rank-sum (Mann–Whitney) test
#'   \item \code{shapiro_test} — p-value from a Shapiro–Wilk normality test applied to the pooled data
#'   \item \code{bartlett_test} — p-value from Bartlett's test for homogeneity of variances
#'   \item \code{levene_test} — p-value from Levene's test for equality of variances
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Identifies all unique pairs of groups in \code{group} (excluding \code{NA})
#'   \item Subsets \code{var} for each group
#'   \item Runs the selected tests
#'   \item Combines results into a single table
#' }
#'
#' Bartlett's test assumes normality and is sensitive to deviations from it.
#' Levene's test is more robust to non-normal data.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   set.seed(123)
#'   values <- rnorm(12)
#'   groups <- rep(c("A", "B", "C"), each = 4)
#'   results <- getTestsForAllCombinations(values, groups)
#'   print(results)
#' }
#' }
#'
#' @seealso
#'  \code{\link[dplyr]{arrange}},
#'  \code{\link[stats]{t.test}},
#'  \code{\link[stats]{wilcox.test}},
#'  \code{\link[stats]{shapiro.test}},
#'  \code{\link[stats]{bartlett.test}},
#'  \code{\link[car]{leveneTest}}
#'
#' @rdname getTestsForAllCombinations
#' @export
#' @importFrom dplyr arrange bind_rows
#' @importFrom stats t.test wilcox.test shapiro.test bartlett.test
#' @importFrom car leveneTest
getTestsForAllCombinations <- function(var, group, sep="-", paired=F){
  comp <- combn(unique(group[!is.na(group)]), 2, simplify = F)
  res <- lapply(comp, function(combo){
    a <- var[group==combo[1]]
    b <- var[group==combo[2]]
    res <- data.frame(groups_compared = paste(combo, sep=sep, collapse=sep),
                      t_pval = t.test(a, b, paired=paired)$p.value,
                      wilcox_pval = wilcox.test(a, b, paired=paired)$p.value,
                      shapiro_test =  shapiro.test(c(a, b))$p.value,
                      bartlett_test = bartlett.test(var[group %in% combo],
                                                    group[group %in% combo])[[3]],
                      levene_test = leveneTest(var[group %in% combo],
                                               factor(group[group %in% combo]))[1, 3]
    )
  }) %>% bind_rows  %>% dplyr::arrange(groups_compared)
  return(res)
}
