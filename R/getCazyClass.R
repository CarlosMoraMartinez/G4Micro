getCazyClass <- function(cazy_tt){
  cazytypes <- c( "GT"="GlycosylTransferases",
                  "GH"="Glycoside Hydrolases",
                  "PL"="Polysaccharide Lyases",
                  "CE"="Carbohydrate Esterases",
                  "AA"="Auxiliary Activities",
                  "CB"="Carbohydrate-Binding Modules")
  split_list <- str_split(cazy_tt, pattern = "\\|") %>% unlist %>% substr(1,2)
  cazy_class <- cazytypes[split_list]
  return(cazy_class)
}
