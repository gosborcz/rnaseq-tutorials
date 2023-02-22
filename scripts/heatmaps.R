# One of the most important part of the RNAseq analysis is visualistaion.
# Take a look at the example code that loead data from another projects and plots a heatmap 
# try to modify it to plot results from all analyses


library(readxl)
library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(gplots)

results <- read_xlsx('/projects/portsmouth-darek-cerebellum/results/edger-results-table.xlsx')

samples <- data.frame(c('sample_BD1', 'sample_BD2', 'sample_BD3', 'sample_BD4', 'sample_BD5', 'sample_BW1', 'sample_BW2', 'sample_BW3', 'sample_BW4',
                        'sample_BW5', 'sample_MD1', 'sample_MD2', 'sample_MD3', 'sample_MD4', 'sample_MD5', 'sample_MW1', 'sample_MW2', 'sample_MW3',
                        'sample_MW4', 'sample_MW5'), c(rep('wt',10), rep('mdx', 10)),
                      c(rep('10days', 5), rep('10weeks', 5), rep('10days', 5), rep('10weeks', 5)),
                      c(rep('WT10D', 5), rep('WT10W', 5), rep('MDX10D', 5), rep('MDX10W', 5)))
colnames(samples) <- c('id', 'genotype', 'age', 'group')
rownames(samples) <- samples$id


#### apply filter for plotting and keep gene names ####

to_plot <- results %>% 
  filter(FDR_genotype < 0.1) %>% 
  select(samples$id)

rownames(to_plot) <- results %>% 
  filter(FDR_genotype < 0.1) %>% 
  select(gene_name) %>% t()


#### actual plotting ####

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- unique(as.character(samples$group[order(samples$id)]))

col.labels <- c(rep("", 2), group.names[1], rep(" ", 4), 
                group.names[2], rep(" ", 4),
                group.names[3], rep(" ", 4),
                group.names[4], rep(" ", 2))


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
    key=TRUE,
    keysize = 0.3,
    Colv = FALSE,
    scale="row",
    colsep = c(5,10,15),
    sepwidth = c(0.3,0.3),
    labRow=rownames(to_plot),
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.7,
    offsetCol = 0,
    margins = c(5, 5)
  )

