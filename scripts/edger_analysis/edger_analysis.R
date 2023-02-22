# installing EdgeR package - run this only once in the R studio

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("edgeR")

# script starts here

#STEP 1 load necessary libraries
require(edgeR) # link to user guide: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
require(tidyverse)

#STEP 4 add group infor th

ctrls <- c("H27", "H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

sample_info <- cbind(
  c(ctrls,ruptured),
  c(rep("ctrl", 4),rep("rup",4)
  )
) %>% 
  as.data.frame() %>% `colnames<-`(c("samples", "groups")) 



#STEP 2 read in the raw data - 
rawdata <- read_tsv("/rnaseq-tutorials/data/with_ref/salmon.merged.gene_counts.tsv")

#STEP 3 load data into edger
df <- DGEList(
  counts = rawdata[,],
  genes = rawdata[,1:2]
)
rownames(df$counts) <- rownames(df$genes) <- df$genes$gene_id


ctrls <- c("H27", "H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

df$samples # see the samples (and the order)

df$samples$group <- c(
  "rup",
  "ctrl",
  "ctrl",
  "rup",
  "rup",
  "ctrl",
  "rup",
  "ctrl"
)

df$samples

keep <- filterByExpr(df)

table(keep)

df <- df[keep, , keep.lib.sizes = FALSE]
df <- calcNormFactors(df)



df <- estimateDisp(df)
et <- exactTest(df, pair = c("ctrl", "rup"))
et_results <- topTags(et, n = Inf)
et_results <- et_results$table %>% tibble()

summary(decideTests(et))

plotMDS(df)


plotMD(et)
abline(h = c(-1, 1), col = "blue")

# Task 2 - add the columns with counts from raw data to the results
# and write the file as tsv

et_results %>%
  cbind(
    rawdata[match(et_results$gene_id, rawdata$gene_id), c(ctrls, ruptured)]
  ) %>%
  write_tsv("/rnaseq-tutorials/data/with_ref/edger_results.csv")

# Task 3 - annotate the results with transcript types from the gtf file from yeasterday.
# How many transcripts of each type are regulated?

# Task 4 - use enrichR (manually or via R) to look for overrepresented pathways in the top results
