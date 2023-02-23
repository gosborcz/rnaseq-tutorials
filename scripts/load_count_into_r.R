require(dplyr)
require(magrittr)

data <- read.delim(
  "/seq/malus_domestica_inra/rnaseq/star_salmon/salmon.merged.gene_counts.tsv",
  header = TRUE
  )

getid <- function(x) {
  x %>% strsplit(":") %>% unlist %>% last
}
gene_id <- data$gene_id %>% lapply(getid) %>% unlist
data <- data[,-c(1,2)] %>% as.matrix
data <- data[,match(sample_sheet$sample, colnames(data))]
