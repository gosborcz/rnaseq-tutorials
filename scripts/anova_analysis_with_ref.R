# a trick to install the preprocessCore library in newer R verions:

BiocManager::install("preprocessCore", configure.args="--disable-threading")

# the script starts here 
require(tidyverse)
require(preprocessCore)

raw_data <- read_tsv("/rnaseq-tutorials/data/with_ref/salmon.merged.gene_tpm.tsv")

ctrls <- c("H27","H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

sample_info <- cbind(
  c(ctrls,ruptured),
  c(rep("ctrl", 4),rep("rup",4)
    )
  ) %>% 
  as.data.frame() %>% `colnames<-`(c("samples", "groups")) 

results <- raw_data %>% 
  select(c("gene_id", "gene_name")) 

tpms <- raw_data[,sample_info$samples]
rownames(tpms) <- raw_data$gene_id

tpms <- data.matrix(tpms)
tpms.norm <-  normalize.quantiles(tpms)
tpms.log <- log2(tpms.norm + 1)

colnames(tpms.log) <- sample_info$samples
results <- cbind(results, tpms.log)

results <- results %>% 
  filter(rowMeans(results[,sample_info$samples]) > 1)

one.way <- function(counts) {
  aov(as.numeric(counts) ~ sample_info$groups) %>%
    summary() %>% unlist() %>% magrittr::extract("Pr(>F)1")
}

apply(results[,sample_info$sample],
      1,
      one.way) %>% 
  p.adjust(method = 'fdr') %>%
  cbind(results, .) -> results

colnames(results)[11] <- "FDR"

results %>% 
  filter(FDR < 0.05) %>% 
  dplyr::select(gene_name) %>% 
  distinct() %>% 
  write.table(row.names = FALSE, quote=FALSE)






