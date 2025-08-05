#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param samples PARAM_DESCRIPTION
#' @param daataxa PARAM_DESCRIPTION
#' @param vstdf PARAM_DESCRIPTION
#' @param net_estimator PARAM_DESCRIPTION
#' @param net_method PARAM_DESCRIPTION
#' @param daatab PARAM_DESCRIPTION
#' @param daatab2 PARAM_DESCRIPTION, Default: NULL
#' @param outdir PARAM_DESCRIPTION, Default: './'
#' @param filter_empty PARAM_DESCRIPTION, Default: FALSE
#' @param filt_quantile PARAM_DESCRIPTION, Default: 0.95
#' @param name PARAM_DESCRIPTION, Default: 'testnet'
#' @param w PARAM_DESCRIPTION, Default: 12
#' @param h PARAM_DESCRIPTION, Default: 12
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select_all}}
#'  \code{\link[assertthat]{assert_that}}
#' @rdname getGraphFromSamples
#' @export 
#' @importFrom dplyr filter mutate select_if
#' @importFrom assertthat assert_that
getGraphFromSamples <- function(phobj, samples, daataxa, vstdf,
                                net_estimator, net_method,
                                daatab, daatab2=NULL,
                                outdir="./",
                                filter_empty=FALSE,
                                filt_quantile=0.95,
                                name="testnet", w=12, h=12){
  # https://github.com/ryanjw/co-occurrence
  # https://rdrr.io/bioc/minet/man/minet.html
  auxf <- function(x){sd(x)>0}
  df2net <- vstdf %>%
    dplyr::filter(gene %in% daataxa) %>%
    dplyr::mutate(gene = process_names(gene)) %>%
    column_to_rownames("gene") %>%
    as.matrix %>% t %>% data.frame %>%
    rownames_to_column("sample") %>%
    dplyr::filter(sample %in% samples) %>%
    column_to_rownames("sample") %>%
    dplyr::select_if(auxf)

  if(net_method == "COR"){
    netres <- cor(df2net, method=net_estimator)
    net_pvals <- map(df2net, \(x)
                     map_vec(df2net, \(y)cor.test(x, y, method = net_estimator)$p.value)) %>%
      data.frame() %>%
      as.matrix
  }else{
    netres <- minet(df2net,
                    method=net_method,
                    estimator=net_estimator,
                    disc="none",
                    nbins=sqrt(NROW(df2net)))
  }




  vert_cols <- ifelse(is.na(daatab$padj), C_NS,
                      ifelse(daatab$padj < opt$pval,
                             ifelse(daatab$log2FoldChangeShrink > 0, C_CASE, C_CTRL), C_NS))
  if(! is.null(daatab2)){
    daatab2 <- daatab2 %>% dplyr::filter(taxon %in% daatab$taxon)
    daatab2 <- daatab2[match(daatab2$taxon, daatab$taxon), ]
    assertthat::assert_that(all(daatab2$taxon == daatab$taxon)) # "Taxa names of daatab1 and daatab2 are not equal."
    vert_cols <- ifelse(is.na(daatab2$padj) | vert_cols != C_NS, vert_cols,
                        ifelse(daatab2$padj < opt$pval,
                               ifelse(daatab2$log2FoldChangeShrink > 0, C_CASE3, C_CTRL3),
                               vert_cols))
  }

  names(vert_cols) <- process_names(daatab$taxon)
  assertthat::assert_that(all(names(df2net) %in% names(vert_cols) ), msg = "Error! names in dfnet and vert colors are different")

  print(table(vert_cols))

  if(filt_quantile > 0 && net_method != "COR"){
    filtthres <- quantile(netres, filt_quantile)
    netres <- ifelse(netres < filtthres, 0, 1)
  }else if(filt_quantile > 0 && net_method == "COR"){
    netres <- ifelse(net_pvals > 1-filt_quantile, 0, 1)
    diag(netres) <- 0
  }
  if(filter_empty){
    keep <- apply(netres, MAR=1, \(x)any(x>0)) %>% which %>% names
    netres <- netres[, keep]
  }

  net <- graph.adjacency(netres, mode="undirected")
  vert_cols <- vert_cols[names(V(net))]
  print(length(V(net)))

  fname <- paste0(outdir, "/",name, "net.pdf")
  pdf(fname, width=12, height=12)
  p1 <- plot(net, vertex.label.size=0.1, vertex.label=NA,
             vertex.size=4, vertex.color=vert_cols)
  dev.off()
  return(list(net=net, cols=vert_cols))
}
