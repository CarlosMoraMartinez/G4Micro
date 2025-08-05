getGraphProps <- function(net_obj){
  net <- net_obj$net
  gprop <- data.frame(
    vertex = names(V(net)),
    closeness = closeness(net),
    betweenness = betweenness(net),
    degree = degree(net),
    color = net_obj$cols
  ) %>%
    dplyr::mutate(closeness = ifelse(is.na(closeness),0,closeness)) %>%
    dplyr::mutate(across(all_of(c("closeness", "betweenness", "degree")), \(x){
      scaled <- (x-min(x))/(max(x)-min(x))
      scaled[is.na(scaled)] <- 0
      return(scaled)
    },.names = "scaled_{.col}")) %>%
    dplyr::mutate(sum_measures = scaled_closeness+scaled_betweenness+scaled_degree) %>%
    dplyr::arrange(desc(sum_measures))
  return(gprop)
}
