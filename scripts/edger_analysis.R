# installing EdgeR package - run this only once

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("edgeR")

# script starts here, remember to add your comments

require(edgeR)
require(tidyverse)

# Task 1 - read the rnaseq-tutorials/data/with_ref/salmon.merged.gene_counts.tsv file into R

rawdata <- # put your answer here

df <- DGEList(
  counts = rawdata[, 3:10],
  genes = rawdata[, 1:2]
)

rownames(df$counts) <- rownames(df$genes) <- df$genes$gene_id

keep <- filterByExpr(df)

table(keep)
df <- df[keep, , keep.lib.sizes = FALSE]
df <- calcNormFactors(df)

plotMDS(df)

ctrls <- c("H27", "H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

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

df <- estimateDisp(df)
et <- exactTest(df, pair = c("ctrl", "rup"))
et_results <- topTags(et, n = Inf)
et_results <- et_results$table %>% tibble()

summary(decideTests(et))
plotMD(et)
abline(h = c(-1, 1), col = "blue")

# Task 2 - add the columns with counts from raw data to the results
# and write the file as tsv

# Task 3 - use enrichR manually to look for overrepresented pathways in the top results
# check upregulated and downregulated fies

# Task 4* - annotate the results with transcript types from the gtf file from yeasterday.
# How many transcripts of each type are regulated?
