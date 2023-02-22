# a trick to install the preprocessCore library in newer R verions
# run ony once; if you get an error in normalize quantiles delete the library restart R and try again

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("preprocessCore", configure.args="--disable-threading")

### the script starts here ###

# STEP 1 load required libraries
require(tidyverse)
require(preprocessCore)

# STEP 2 read in data - remember for this type of analysis use tmp files
# remember to change the paths to your own
raw_data <- read_tsv(
  "/rnaseq-tutorials/data/with_ref/salmon.merged.gene_tpm.tsv"
  )

# STEP 3 make two vectors (manually) with control and experimental samples
ctrls <- c("H27","H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

# STEP 4 make a table with sample info 
sample_info <- cbind(
  c(ctrls,ruptured), # combine the sample names from vectors defined above
  c(rep("ctrl", 4),rep("rup",4) # make a new column with the names of the group
    )
  ) %>% 
  as.data.frame() %>% `colnames<-`(c("samples", "groups")) 


# STEP 5 make a new table that will host results, start with the gene_id and gene_name columns
results <- raw_data %>% 
  select(c("gene_id", "gene_name")) 

# STEP 6 subset the raw data to only have counts (numbers) and set its rownames to gene_id
tpms <- raw_data[,sample_info$samples]
rownames(tpms) <- raw_data$gene_id

# STEP 7 make a matrix of numbers from the tpms, normalize it
# and give it sample names as column names 
tpms <- data.matrix(tpms)
tpms.norm <-  normalize.quantiles(tpms)
tpms.log <- log2(tpms.norm + 1)
colnames(tpms.log) <- sample_info$samples

# STEP 8 join the results table created above with the log-transformed counts
results <- cbind(results, tpms.log)

# STEP 9 filter out lowly expressed genes
results <- results %>% 
  filter(rowMeans(results[,sample_info$samples]) > 1)

# STEP 10 make a one-way anova fuction that extracts the p-value
one.way <- function(counts) {
  aov(as.numeric(counts) ~ sample_info$groups) %>%
    summary() %>% unlist() %>% magrittr::extract("Pr(>F)1")
}

# STEP 11 apply the function above to all rows
by_rows <- 1 #number 1 represents rows

apply(results[,sample_info$sample], # this applys the one-way function
      by_rows,
      one.way) %>% 
  p.adjust(method = 'fdr') %>% # this adjusts the p-value
  cbind(results, .) -> results # this binds the column with p-values to results

# STEP 12 rename the column to FDR
colnames(results)[11] <- "FDR" 

# STEP 13 the function below calculates group means
compute_means <- function(df, sample_info, sample_names_re, factors) {
  sample_info_grouped <-
    sample_info %>%
    group_by(across(all_of(factors)))
  groups <-
    sample_info_grouped %>%
    group_split()
  group_keys <-
    sample_info_grouped %>%
    group_keys() %>% 
    unite('name', everything(), remove = FALSE)
  
  for (tb_idx in 1:length(groups)) {
    tb <- groups[[tb_idx]]
    tb_key <- group_keys[tb_idx, 'name'] %>% unlist()
    df <-
      df %>%
      mutate(
        "{tb_key}_mean" := apply(
          select(., matches(sample_names_re)),
          1,
          function(x, names) {mean(x[names])},
          names = tb$id
        )
      )
  }
  df
}

# STEP 14 run the function to compute means
# to run the function above we need a column in sample_info named id
sample_info$id <- sample_info$samples

compute_means(
  results[,sample_info$samples],
  sample_info,
  'H', # this is a pattern that finds sample names in results
  factors = c('groups') # this is the name of the groups column in sample info
) -> group_means # this creates a new table group_means where last columns are mean log2FC for each group


# STEP 15 add the log2 fold change based on the group_means to the results
results$log2FC <- group_means$rup_mean-group_means$ctrl_mean # negative values indicates downregulation and positive - upregulation

# STEP 16 print a list of top genes
results %>% 
  filter(FDR < 0.05) %>% # FDR < 5% - can be adjusted
  filter(log2FC > 1) %>% # direction UP- can be adjusted
  dplyr::select(gene_name) %>%
  write.table(row.names = FALSE, quote=FALSE)

# STEP 17 save the results
results %>% 
  write_tsv("/rnaseq-tutorials/data/with_ref/anova_results.csv")






