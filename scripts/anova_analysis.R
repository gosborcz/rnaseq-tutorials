# a trick to install the preprocessCore library in newer R verions:

BiocManager::install("preprocessCore", configure.args="--disable-threading")

# the script starts here 
require(plyr)
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

results <- dfs %>% 
  dplyr::select(c("Name", "gene_names")) 

# Task 2 - find out what does the dplyr:: before the function do

tpms <- dfs[,sample_info$samples]
rownames(tpms) <- dfs$Name

tpms <- data.matrix(tpms)
tpms.norm <-  normalize.quantiles(tpms)

# Task 3 - Take a look at the boxplots of TPMs before and after normalization

tpms.log <- log2(tpms.norm + 1)

colnames(tpms.log) <- sample_info$samples
results <- cbind(results, tpms.log)

# Task 4 - replace rows with tpms.log that have row means < 1

results <- results %>% 
  filter(rowMeans(results[,sample_info$samples]) > 1)

results <- results %>% 
  filter(!is.na(results$gene_names)) %>% 
  pivot_longer(sample_info$samples,
                         names_to = "sample",
                         values_to = "log2tpm") %>%
  group_by(gene_names,sample) %>%
  dplyr::summarize(mean = mean(log2tpm)) %>%
  pivot_wider(id_cols = gene_names,
              names_from = sample,
              values_from = mean)

one.way <- function(counts) {
  aov(as.numeric(counts) ~ sample_info$groups) %>%
    summary() %>% unlist() %>% extract("Pr(>F)1")
}

# Task 5 apply this function to the log2-transformed tpms

apply(results[,sample_info$sample],
      1,
      one.way) %>% 
  #p.adjust(method = 'fdr') %>%
  cbind(results, .) -> results

colnames(results)[10] <- "FDR"

results %>% 
  filter(FDR < 0.1) %>% 
  dplyr::select(gene_names) %>% 
  distinct() %>% 
  write.table(row.names = FALSE, quote=FALSE)

# Task 6 calculate log2 fold-change (log2FC) for each detected transcript

# Task 7 = perfom an analysis with one-way anova of the results 
# from the alignment to the reference genome in a separate script
# are the results different from Edger?

# Task 8* - downolad data from this publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7691544/
# raw counts are here: http://149.156.177.112/projects/ifpan-kinga-dieta/cuffnorm/genes.count_table
# sample info here: http://149.156.177.112/projects/ifpan-kinga-dieta/cuffnorm/samples.table
# and here: http://149.156.177.112/projects/ifpan-kinga-dieta/analysis/samples-info.csv
#  and run your own two-factor analysis in a new file






