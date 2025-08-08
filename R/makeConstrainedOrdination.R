#'  @title Make a Constrained Ordination Plot (CAP) from Phyloseq Object
#' @description Performs a constrained ordination (Canonical Analysis of Principal coordinates - CAP)
#'   on a phyloseq object using specified metadata variables to constrain the ordination.
#'   It plots the ordination with environmental vectors (arrows) representing constraining variables.
#'
#' @param phobj A phyloseq object containing microbiome data with sample metadata.
#' @param variables Character vector of metadata variable names to constrain the ordination by.
#'   Default is \code{c("Condition")}.
#' @param dist_type Character string specifying the distance metric to use for sample dissimilarity.
#'   Default is \code{"bray"}.
#' @param color_by Character string specifying the metadata variable to use for coloring the points
#'   in the ordination plot. If empty, defaults to the first variable in \code{variables}.
#'
#' @return A ggplot2 object showing the constrained ordination plot with samples colored by \code{color_by},
#'   and environmental vectors (arrows) for each constraining variable.
#'
#' @details
#' This function removes samples with missing values in any of the constraining variables,
#' calculates a distance matrix using the specified distance metric, performs CAP ordination,
#' and visualizes the results including arrows for environmental variables.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(phyloseq)
#'   # Assuming ps is a phyloseq object with "Condition" metadata
#'   makeConstrainedOrdination(ps, variables = c("Condition", "Timepoint"), dist_type = "bray", color_by = "Condition")
#' }
#' }
#'
#' @seealso
#'  \code{\link[phyloseq]{distance}}
#'  \code{\link[phyloseq]{ordinate}}
#'  \code{\link[phyloseq]{plot_ordination}}
#'  \code{\link[vegan]{scores}}
#'
#' @rdname makeConstrainedOrdination
#' @export
#' @importFrom phyloseq distance ordinate plot_ordination prune_samples sample_data
#' @importFrom vegan scores
#' @importFrom rlang sym
#' @importFrom ggplot2 aes geom_point geom_segment geom_text
#' @importFrom ggplot2 arrow unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom assertthat assert_that
makeConstrainedOrdination<-function(phobj, variables = c("Condition"),  dist_type="bray", color_by=""){
  assertthat::assert_that(length(variables) > 1)
  phobj_not_na <- phobj
  metad <- sample_data(phobj_not_na) %>% data.frame
  na_samples_2Remove <- metad %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(variables), is.na)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(any_na = any(dplyr::c_across(all_of(variables)))) %>%
    dplyr::ungroup() %>% dplyr::filter(!any_na) %>%
    dplyr::pull(sampleID) # samples to preserve
  phobj_not_na <- prune_samples(na_samples_2Remove, phobj_not_na)
  metad <- sample_data(phobj_not_na) %>% data.frame

  form <- paste0("~ ", paste(variables, sep=" + ", collapse = " + ")) %>% as.formula

  bray_not_na <- phyloseq::distance(physeq = phobj_not_na, method = dist_type)

  # CAP ordinate
  cap_ord <- phyloseq::ordinate(physeq = phobj_not_na, method = "CAP",
                                distance = bray_not_na, formula = form)

  # CAP plot
  if(color_by == ""){
    color_by <- variables[1]
  }
  cap_plot <- phyloseq::plot_ordination(physeq = phobj_not_na, ordination = cap_ord,
                                        color = color_by, axes = c(1,2)) +
    aes(col = !!sym(color_by)) +
    geom_point( size = 4)
  #geom_point(colour = "grey90", size = 1.5, type=1)
  #scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a",
  #                              "#1919ff", "darkorchid3", "magenta"))

  # if(! is.numeric(metad[, color_by])){
  #   cap_plot <- cap_plot + scale_color_lancet()
  # }
  # Now add the environmental variables as arrows
  arrowmat <- vegan::scores(cap_ord, display = "bp")

  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  text_colors <- rep("gray19", length(variables))
  if(color_by %in% variables){
    text_colors[color_by == variables] <- "green4"
  }

  # Define the arrow aesthetic mapping
  arrow_map <- aes(xend = CAP1,
                   yend = CAP2,
                   x = 0,
                   y = 0,
                   shape = NULL,
                   color = NULL,
                   label = labels)

  label_map <- aes(x = 1.3 * CAP1,
                   y = 1.3 * CAP2,
                   shape = NULL,
                   color = NULL,
                   label = labels)

  arrowhead = arrow(length = unit(0.02, "npc"))

  # Make a new graphic
  g1 <- cap_plot +
    geom_text_repel( aes(label=sampleID))+
    geom_segment(
      mapping = arrow_map,
      size = 1,
      data = arrowdf,
      color = text_colors, #"gray",
      arrow = arrowhead
    ) +
    geom_text(
      mapping = label_map,
      size = 6,
      color = text_colors,
      data = arrowdf,
      show.legend = FALSE
    ) + mytheme

  return(g1)
}
