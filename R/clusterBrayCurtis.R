#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param variable PARAM_DESCRIPTION, Default: 'Condition'
#' @param renamestr PARAM_DESCRIPTION, Default: ''
#' @param dist_type PARAM_DESCRIPTION, Default: 'bray'
#' @param clust_method PARAM_DESCRIPTION, Default: 'ward.D2'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[phyloseq]{transform_sample_counts}}, \code{\link[phyloseq]{otu_table}}, \code{\link[phyloseq]{distance}}, \code{\link[phyloseq]{sample_data}}
#' @rdname clusterBrayCurtis
#' @export 
#' @importFrom phyloseq transform_sample_counts otu_table distance sample_data
clusterBrayCurtis <- function(phobj, variable="Condition", renamestr="",
                              dist_type="bray", clust_method="ward.D2"){
  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})
  ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
  ps_rel_otu <- t(ps_rel_otu)
  rownames(ps_rel_otu) <- gsub("^X", renamestr, rownames(ps_rel_otu))
  #bc_dist <- vegan::vegdist(ps_rel_otu, method = dist_type)
  bc_dist <- phyloseq::distance(physeq = phobj, method = dist_type)

  #Save as dendrogram
  ward <- as.dendrogram(hclust(bc_dist, method = clust_method))
  #Provide color codes
  meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
  levv <- unique(meta[, variable])
  colorCode <- c( "firebrick3", "dodgerblue3")
  names(colorCode) <- levv
  labels_colors(ward) <- colorCode[meta[, variable]][order.dendrogram(ward)]

  hmplot <- plot_heatmap(ps_rel_abund, sample.label=variable)
  return(list("ward"=ward, "hmplot"=hmplot))
}
