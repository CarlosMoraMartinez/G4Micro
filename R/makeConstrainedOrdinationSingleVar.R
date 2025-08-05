#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param variables PARAM_DESCRIPTION, Default: c("Condition")
#' @param dist_type PARAM_DESCRIPTION, Default: 'bray'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[phyloseq]{distance}}
#'  \code{\link[vegan]{scores}}
#' @rdname makeConstrainedOrdinationSingleVar
#' @export 
#' @importFrom phyloseq distance
#' @importFrom vegan scores
makeConstrainedOrdinationSingleVar<-function(phobj, variables = c("Condition"),  dist_type="bray"){
  phobj_not_na <- phobj
  metad <- data.frame(sample_data(phobj_not_na))
  na_samples_2Remove <- metad$sampleID[! is.na(metad[, variables[1]])]
  phobj_not_na <- prune_samples(na_samples_2Remove, phobj_not_na)
  #phobj_not_na <- subset_samples(phobj_not_na, ! sampleID %in% na_samples_2Remove)

  form <- paste0("~ ", paste(variables, sep=" + ", collapse = " + ")) %>% as.formula
  bray_not_na <- phyloseq::distance(physeq = phobj_not_na, method = dist_type)

  # CAP ordinate
  cap_ord <<- ordinate(
    physeq = phobj_not_na,
    method = "CAP",
    distance = bray_not_na,
    formula = form
  )

  # CAP plot
  cap_plot <- plot_ordination(
    physeq = phobj_not_na,
    ordination = cap_ord,
    color = variables[1],
    axes = c(1,2)
  ) +
    aes_string(col = variables[1]) +
    geom_point( size = 4)
  #geom_point(colour = "grey90", size = 1.5, type=1)

  # if(! is.numeric(metad[, variables[1]])){
  #   cap_plot <- cap_plot + scale_color_lancet()
  # }

  #scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a",
  #                              "#1919ff", "darkorchid3", "magenta"))


  # Now add the environmental variables as arrows
  arrowmat <- vegan::scores(cap_ord, display = "bp")

  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

  # Define the arrow aesthetic mapping
  arrow_map <- aes(xend = CAP1,
                   yend = 0,
                   x = 0,
                   y = 0,
                   shape = NULL,
                   color = NULL,
                   label = labels)

  label_map <- aes(x = 1.3 * CAP1,
                   y = 1.3 * 0,
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
      color = "gray",
      arrow = arrowhead
    ) +
    geom_text(
      mapping = label_map,
      size = 6,
      data = arrowdf,
      show.legend = FALSE
    ) + mytheme

  return(g1)
}
