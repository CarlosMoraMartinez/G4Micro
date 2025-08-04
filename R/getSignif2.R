getSignif2 <- function(vector, limits=c(0.05, 0.01, 0.001),
                       legends = c("", "*", "**", "***")){
  x <- ifelse(vector < limits[3], legends[4],
              ifelse(vector < limits[2], legends[3],
                     ifelse(vector< limits[1], legends[2], legends[1])
              )
  )
  return(x)
}
