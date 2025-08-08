

options(ggplot2.continuous.fill="viridis")
options(ggplot2.continuous.colour="viridis")

.onLoad <- function(libname, pkgname) {
  options(ggplot2.discrete.fill = c("#A1C6EA","#FD8B2F", "#00AA5A", "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))
  options(ggplot2.discrete.colour = c("#A1C6EA","#FD8B2F","#00AA5A",   "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))
  opt <- opt_default
  restaurar <- restauraropt_mk(opt)
  packageStartupMessage("G4Micro package loaded!")
}
#mycols <- grDevices::colorRampPalette( wesanderson::wes_palette("Royal1"))(5)

# options(ggplot2.discrete.fill = mycols)
# options(ggplot2.discrete.colour = mycols)
