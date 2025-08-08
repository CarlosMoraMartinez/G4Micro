
# Core
SEED <- 123
MULTI_PAGE_PDFS = TRUE

# Colors
C_CASE = "#FD8B2F"
C_CASE2 = "tomato"
C_CASE_LINK = "#fBd895"
C_CTRL = "#A1C6EA"
C_CTRL2 = "steelblue2"
C_CTRL_LINK ="#DAE3E5"
C_CTRL_LINK2 ="#B8C1D4"
C_WHITE= "#DDDDDD"
C_NS =  "#A5ABBD"
C_OTHER = "gray30"

C_CASE3 = "#00AA5A"
C_CTRL3 = "#8E7BFF"

pal1 <- wesanderson::wes_palette("Darjeeling2", 2)
pal <- wesanderson::wes_palette("AsteroidCity2")
C_POS_EFF = pal1[2]
C_NEG_EFF = pal[4]

mycols <- grDevices::colorRampPalette( wesanderson::wes_palette("Royal1"))(5)

# Themes
mytheme <-  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5,
                                  colour = "black", face = "bold")) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.text.y = element_text(size = 12,
                                   colour = "black", angle = 0, face = "bold")) +
  theme(strip.text.y = element_text(size = 12,
                                    colour = "black", angle = 0, face = "bold")) +
  theme(strip.text.x = element_text(size = 12,
                                    colour = "black", angle = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 11,
                                   colour = "black", angle = 0,
                                   face = "bold"))+
  theme(axis.title.x = element_text(vjust = 1, hjust = 0.5,
                                    size = 12, colour = "black",
                                    angle = 0, face = "bold")) +
  theme(axis.title.y= element_text(vjust = 1, hjust = 0.5,
                                   size = 12, colour = "black",
                                   angle = 90, face = "bold"))


mystyle <- theme_classic() +
  theme(axis.text = element_text(face="bold"),
        axis.title = element_text(face="bold"))


thin_barplot_lines <- theme(panel.grid.major.y = element_line(color = "lightgray",
                                                              size = 0.05,
                                                              linetype = 2))


#' @export
opt_default <- list(out ="./",
            minfreq = 0.05,
            mincountspersample = 0,
            mincount= 1,
            minsampleswithcount = 0,
            raref_quant = 0.15,
            fc=1,
            pval=0.05,
            ptype="adjusted",
            fctype="shrunk",
            num_genes_default=5 # meaning genes or taxa, depending on the context
)
