#' @title makeConstrainedOrdinationSingleVar
#' @description Performs a constrained ordination (CAP) on a phyloseq object using specified variables.
#' Uses the \code{\link[phyloseq]{ordinate}} method with the CAP method.
#' @param phobj A phyloseq object containing microbiome data and sample metadata.
#' @param variables A character vector of variables to include in the ordination formula. Default is \code{c("Condition")}.
#' @param dist_type Distance metric to use for ordination. Default is \code{"bray"}.
#' @return A \code{ggplot2} object representing the constrained ordination plot with environmental variables as arrows.
#' @details
#' The function removes samples with NA in the first variable specified, calculates the distance matrix,
#' performs CAP ordination constrained by the variables, and plots the ordination with points colored by the first variable.
#' Environmental variables are added as arrows with labels.
#' Note: This function is redundant with \code{makeConstrainedOrdination}, but this one allows to plot CAP with only 1 variable, whereas the other one needs more than 1.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example assuming `physeq` is a phyloseq object with metadata variable "Condition"
#'   makeConstrainedOrdinationSingleVar(physeq, variables = c("Condition"))
#' }
#' }
#' @seealso
#' \code{\link[phyloseq]{distance}}, \code{\link[vegan]{scores}}, \code{\link[phyloseq]{ordinate}}
#' @rdname makeConstrainedOrdinationSingleVar
#' @export
#' @importFrom phyloseq distance
#' @importFrom vegan scores
#' @importFrom ggrepel geom_text_repel
makeConstrainedOrdinationSingleVar<-function(phobj, variables = c("Condition"),  dist_type="bray"){
  phobj_not_na <- phobj
  metad <- data.frame(sample_data(phobj_not_na))
  na_samples_2keep <- metad$sampleID[! is.na(metad[, variables[1]])]
  phobj_not_na <- prune_samples(na_samples_2keep, phobj_not_na)

  form <- paste0("~ ", paste(variables, sep=" + ", collapse = " + ")) %>% as.formula
  bray_not_na <- phyloseq::distance(physeq = phobj_not_na, method = dist_type)

  # CAP ordinate
  # <<- removed
  cap_ord <- ordinate(
    physeq = phobj_not_na,
    method = "CAP",
    distance = bray_not_na,
    formula = form
  )

  # Now add the environmental variables as arrows
  arrowmat <- vegan::scores(cap_ord, display = "bp")

  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  name2 <- colnames(arrowmat)[2]
  # CAP plot
  cap_plot <- plot_ordination(
    physeq = phobj_not_na,
    ordination = cap_ord,
    color = variables[1],
    axes = c(1,2)
  ) +
    aes(col = !!sym(variables[1])) +
    geom_point( size = 4)

  # Define the arrow aesthetic mapping
  arrow_map <- aes(xend = CAP1,
                   yend = !!sym(name2),
                   x = 0,
                   y = 0,
                   shape = NULL,
                   color = NULL,
                   label = labels)

  label_map <- aes(x = 1.3 * CAP1,
                   y = 1.3 * !!sym(name2),
                   shape = NULL,
                   color = NULL,
                   label = labels)

  arrowhead = arrow(length = unit(0.02, "npc"))

  # Make a new graphic
  g1 <- cap_plot +
    ggrepel::geom_text_repel( aes(label=sampleID))+
    geom_segment(
      mapping = arrow_map,
      linewidth = 1,
      data = arrowdf,
      color = "gray",
      arrow = arrowhead
    ) +
    geom_text(
      mapping = label_map,
      linewidth = 6,
      data = arrowdf,
      show.legend = FALSE
    ) + mytheme

  return(g1)
}
