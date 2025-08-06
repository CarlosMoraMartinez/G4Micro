#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param variables PARAM_DESCRIPTION, Default: c("Condition")
#' @param dist_type PARAM_DESCRIPTION, Default: 'bray'
#' @param color_by PARAM_DESCRIPTION, Default: ''
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
#' @rdname makeConstrainedOrdination
#' @export
#' @importFrom phyloseq distance
#' @importFrom vegan scores
makeConstrainedOrdination<-function(phobj, variables = c("Condition"),  dist_type="bray", color_by=""){
  phobj_not_na <- phobj
  metad <- sample_data(phobj_not_na) %>% data.frame
  for(v in variables){
    na_samples_2Remove <- metad$sampleID[! is.na(metad[, v]) ]
    phobj_not_na <- prune_samples(na_samples_2Remove, phobj_not_na)
    #phobj_not_na <- subset_samples(phobj_not_na, sampleID %in% na_samples_2Remove)
    metad <- sample_data(phobj_not_na) %>% data.frame
  }
  form <- paste0("~ ", paste(variables, sep=" + ", collapse = " + ")) %>% as.formula
  bray_not_na <- phyloseq::distance(physeq = phobj_not_na, method = dist_type)

  # CAP ordinate
  cap_ord <- ordinate(
    physeq = phobj_not_na,
    method = "CAP",
    distance = bray_not_na,
    formula = form
  )

  # CAP plot
  if(color_by == ""){
    color_by <- variables[1]
  }
  cap_plot <- plot_ordination(
    physeq = phobj_not_na,
    ordination = cap_ord,
    color = color_by,
    axes = c(1,2)
  ) +
    aes_string(col = color_by) +
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
