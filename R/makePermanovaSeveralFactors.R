#' @title Perform PERMANOVA Tests for Multiple Factor Models on Phyloseq Data
#' @description
#' Computes PERMANOVA (adonis2) tests using the specified distance metric for multiple models defined by combinations of metadata variables on a phyloseq object.
#'
#' @param phobj A \code{phyloseq} object containing microbiome data and sample metadata.
#' @param dist_method Character string specifying the distance metric for beta diversity calculation. Default is "bray".
#' @param seed Integer random seed for reproducibility. Default is 123.
#' @param modlist A list of character vectors, each vector specifying the metadata variables to include in a PERMANOVA model. Example: list(c("Condition"), c("Condition", "Batch")).
#' @param outname File path to save the resulting tibble with PERMANOVA results as an \code{RData} file. Default is "permanovas_mult.RData".
#'
#' @return A tibble containing the formulas tested, the fitted PERMANOVA models (\code{adonis2} output as a list column), and the proportion of variance explained by the model (1 - residual R2).
#'
#' @details
#' For each set of variables in \code{modlist}, a PERMANOVA model is fit with the
#' formula \code{distance ~ variables} using \code{vegan::adonis2}.
#' The function saves the results in an \code{RData} file and returns the tibble for immediate use.
#'
#' Based on: \url{From https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova}
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(phyloseq)
#'   data(GlobalPatterns)
#'   sample_data(GlobalPatterns)$Condition <- sample(c("A", "B"), nsamples(GlobalPatterns), replace = TRUE)
#'   sample_data(GlobalPatterns)$Batch <- sample(c("X", "Y"), nsamples(GlobalPatterns), replace = TRUE)
#'
#'   results <- makePermanovaSeveralFactors(
#'     phobj = GlobalPatterns,
#'     dist_method = "bray",
#'     modlist = list(c("Condition"), c("Condition", "Batch")),
#'     outname = "permanova_results.RData"
#'   )
#'   print(results)
#' }
#' }
#'
#' @seealso
#' \code{\link[phyloseq]{distance}}, \code{\link[vegan]{adonis2}}
#'
#' @rdname makePermanovaSeveralFactors
#' @export
#' @importFrom phyloseq distance sample_data
#' @importFrom vegan adonis2
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
makePermanovaSeveralFactors <- function(phobj,
                                        dist_method = "bray",
                                        seed = 123,
                                        modlist = c(),
                                        outname = "permanovas_mult.RData"){

  meta <- sample_data(phobj)
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))

  res <- tibble()
  for(vars in modlist){
    set.seed(seed)
    form <- paste0("braydist ~ ", paste(vars, sep= " + ", collapse = " + "))
    mod1 <- adonis2(as.formula(form), data = sampledf, na.action=na.exclude)
    aux <- tibble(formula = form, model = list(mod1), explained=1 - mod1["Residual", "R2"])
    res <-rbind(res, aux)
  }
  save(res, file=outname)
  return(res)

}
