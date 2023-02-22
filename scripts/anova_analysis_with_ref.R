# a trick to install the preprocessCore library in newer R verions:

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("preprocessCore", configure.args="--disable-threading")

# the script starts here 
require(tidyverse)
require(preprocessCore)

#task 1 read the data/with_ref/salmon.merged.gene_tpm.tsv file into R 

ctrls <- c("H27","H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

sample_info <- cbind(
  c(ctrls,ruptured),
  c(rep("ctrl", 4),rep("rup",4)
    )
  ) %>% 
  as.data.frame() %>% `colnames<-`(c("samples", "groups")) 

# task 2 - create a new table named results 
# that will have only the gene_id and gene_name columns 

results <- raw_data %>% 
  select(c("gene_id", "gene_name")) 

tpms <- raw_data[,sample_info$samples]
rownames(tpms) <- raw_data$gene_id

tpms <- data.matrix(tpms)
tpms.norm <-  normalize.quantiles(tpms)

# task 3 - create a new matrix named tpms.log
# that will have tpms.norm but normalised according to the log2(x+1) formula


colnames(tpms.log) <- sample_info$samples

#task 4 - join the results table with the tmps.log matrix
#task 5 - filter the restuls table to only have genes with row means > 1

one.way <- function(counts) {
  aov(as.numeric(counts) ~ sample_info$groups) %>%
    summary() %>% unlist() %>% magrittr::extract("Pr(>F)1")
}

# task 6 - apply the one.way function to the valuse in the results table 

# input the code here put a %>% at the end
  p.adjust(method = 'fdr') %>%
  cbind(results, .) -> results

colnames(results)[11] <- "FDR"

results %>% 
  filter(FDR < 0.05) %>% 
  dplyr::select(gene_name) %>% 
  distinct() %>% 
  write.table(row.names = FALSE, quote=FALSE)

results %>% 
  write_tsv("/rnaseq-tutorials/data/with_ref/edger_results.csv")






