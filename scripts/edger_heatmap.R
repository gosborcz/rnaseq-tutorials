library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(gplots)

# Read results from edger analysis
results <- read_tsv("/rnaseq-tutorials/data/with_ref/edger_results.csv")

ctrls <- c("H27", "H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

sample_info <- cbind(
  c(ctrls,ruptured),
  c(rep("ctrl", 4),rep("rup",4)
  )
) %>% 
  as.data.frame() %>% `colnames<-`(c("samples", "groups")) 


# apply filters for plotting and keep gene names ##

to_plot <- results %>% 
  filter(FDR < 0.0001) %>% 
  select(sample_info$samples)

rownames(to_plot) <- results %>% 
  filter(FDR < 0.0001) %>% 
  select(gene_name) %>% t()

#### actual plotting ####

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- unique(as.character(sample_info$groups))

col.labels <- c(rep("", 1), group.names[1], rep(" ", 4),  group.names[2])


cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}


to_plot %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to_plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    main="",
    Colv = FALSE,
    scale="row",
    colsep = c(4),
    sepwidth = c(0.1,0.1),
    labRow=rownames(to_plot),
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.7,
    offsetCol = 0
  )
