nameProcesses<- function(tab, namedtab){
  if(is.null(namedtab)){return(tab)}
  tab$process_name <- namedtab$long[match(rownames(tab), namedtab$short)]
  return(tab)
}
