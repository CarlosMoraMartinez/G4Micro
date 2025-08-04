getSignif <- function(vector){
  x <- ifelse(vector < 0.001, "***",
              ifelse(vector < 0.01, "**",
                     ifelse(vector< 0.05, "*", "NS")
              )
  )
  return(x)
}
