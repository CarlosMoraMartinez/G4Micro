#' @title Run PERMANOVA on All Metadata Variables in a Phyloseq Object
#'
#' @description
#' Runs PERMANOVA (\code{\link[vegan]{adonis2}}) on all metadata variables in a
#' \code{\link[phyloseq]{phyloseq}} object, excluding those specified in
#' \code{exclude_vars}. For each variable, the function:
#' \enumerate{
#'   \item Filters samples with non-missing values.
#'   \item Computes the selected distance matrix (e.g., Bray–Curtis, Jaccard).
#'   \item Performs PERMANOVA to assess group differences.
#'   \item Optionally tests for differences in multivariate dispersion
#'         (\code{\link[vegan]{betadisper}} and \code{\link[vegan]{permutest}}).
#'   \item Optionally tests the significance of constrained ordination
#'         (\code{\link[vegan]{capscale}} with \code{\link[vegan]{anova.cca}}).
#' }
#' Results are adjusted for multiple testing (FDR) and written to a TSV file.
#'
#' @param phobj A \code{\link[phyloseq]{phyloseq}} object containing the OTU/ASV
#'   table and sample metadata.
#' @param dist_method Character string specifying the dissimilarity index to
#'   use (e.g., \code{"bray"}, \code{"jaccard"}). Passed to
#'   \code{\link[phyloseq]{distance}}. Default is \code{"bray"}.
#' @param seed Integer random seed for reproducibility of PERMANOVA results.
#'   Default is \code{123}.
#' @param exclude_vars Character vector of metadata variable names to exclude
#'   from PERMANOVA testing. Default is \code{c("sampleID")}.
#' @param outname Character string with the file name to write results table
#'   to (TSV format). Default is \code{"permanovas.tsv"}.
#' @param disp_permutations Integer number of permutations to use in
#'   PERMANOVA, dispersion tests, and ordination significance tests.
#'   Default is \code{1000}.
#'
#' @return A \code{data.frame} with PERMANOVA results for each tested metadata
#'   variable, including:
#'   \itemize{
#'     \item \code{P} – raw PERMANOVA p-value
#'     \item \code{padj} – adjusted p-value (BH FDR)
#'     \item Optional columns for dispersion and capscale ANOVA p-values and their adjustments
#'   }
#'
#' @details
#' This function is particularly useful for exploratory screening of
#' associations between microbial community composition and metadata variables.
#' It is based on the workflow described at:
#' \url{https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova}
#'
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   makePermanova(
#'     phobj,
#'     dist_method = "bray",
#'     seed = 123,
#'     exclude_vars = c("sampleID", "patient_name"),
#'     outname = "my_output_name.tsv"
#'   )
#' }
#' }
#' @seealso
#' \code{\link[phyloseq]{distance}},
#' \code{\link[phyloseq]{prune_samples}},
#' \code{\link[phyloseq]{sample_data}},
#' \code{\link[vegan]{adonis2}},
#' \code{\link[vegan]{betadisper}},
#' \code{\link[vegan]{permutest}},
#' \code{\link[vegan]{capscale}},
#' \code{\link[vegan]{anova.cca}},
#' \code{\link[dplyr]{arrange}},
#' \code{\link[readr]{write_tsv}}
#'
#' @rdname makePermanova
#' @export
#'
#' @importFrom vegan permutest adonis2 betadisper capscale anova.cca
#' @importFrom phyloseq distance prune_samples sample_data
#' @importFrom dplyr arrange
#' @importFrom stringi stri_trans_general
#' @importFrom stats as.formula p.adjust na.exclude
#' @importFrom readr write_tsv
makePermanova <- function(phobj, dist_method = "bray", seed = 123,
                          exclude_vars = c("sampleID"), outname = "permanovas.tsv", disp_permutations=1000){

  meta <- sample_data(phobj)
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))
  exclude_vars <- stri_trans_general(str = exclude_vars, id = "Latin-ASCII") %>%
    gsub(" ", "_", .)
  names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>%
    gsub(" ", "_", .)

  vars2test <- names(sampledf)[! names(sampledf) %in% exclude_vars]

  res <- data.frame()
  for(var in vars2test){
    set.seed(seed)
    form <- paste0("braydist ~ ", var) %>% as.formula

    if(length(unique(sampledf[!is.na(sampledf[, var]) , var])) > 1){
      # Adonis test
      #cat(var, "\n")
      # get original data frame
      sampledf <- data.frame(sample_data(phobj))
      names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>%
        gsub(" ", "_", .)

      s2use <- sampledf$sampleID[!is.na(sampledf[, var]) ]
      phobj_filt <- phyloseq::prune_samples(s2use, phobj)

      braydist <- phyloseq::distance(phobj_filt, method = dist_method)
      sampledf <- data.frame(sample_data(phobj_filt))
      names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>%
        gsub(" ", "_", .)

      mod1 <- adonis2(form, data = sampledf, na.action=na.exclude, permutations = disp_permutations)

      perm_disp <- tryCatch({
        bd <- betadisper(braydist, sampledf[, var])
        permutest(bd, permutations = disp_permutations)
      },error=\(x) NULL )

      ancap <- tryCatch({
        cap <- capscale(form, data = sampledf)
        ancap <- anova.cca(cap, permutations = disp_permutations)
      },error=\(x) NULL )

      res <-rbind(res, adonis2table(mod1, perm_disp, ancap, var))
    }
  }
  res <- res %>%
    dplyr::arrange(P)
  res$padj <- p.adjust(res$P, method="BH")
  if("perm_disp_P" %in% names(res)){
    res$perm_disp_Padj <- p.adjust(res$perm_disp_P, method="BH")
  }
  if("capscaleanova_P" %in% names(res)){
    res$capscaleanova_Padj <- p.adjust(res$capscaleanova_P, method="BH")
  }
  write_tsv(res, file=outname)
  return(res)

}
