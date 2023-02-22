# installing EdgeR package - run this only once in the R studio

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("edgeR")

# actual script starts here

# STEP 1 load necessary libraries
require(edgeR) # link to user guide: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
require(tidyverse)

# STEP 2 create vectors with group names and make a table with sample info
ctrls <- c("H27", "H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

sample_info <- cbind(
  c(ctrls,ruptured),
  c(rep("ctrl", 4),rep("rup",4)
  )
) %>% 
  as.data.frame() %>% `colnames<-`(c("samples", "groups")) 


# STEP 3 read in the raw data - edgeR need the file with RAW counts
rawdata <- read_tsv("/rnaseq-tutorials/data/with_ref/salmon.merged.gene_counts.tsv")

# STEP 4 load data into edger (the DGEList is an EdgeR function to load data)
df <- DGEList(
  counts = rawdata[,sample_info$samples], #provide columns with counts
  genes = rawdata[,c("gene_id","gene_name")], #provide columns with genes
  group = sample_info$groups #provide group info in the same order as columns with counts
)
rownames(df$counts) <- rownames(df$genes) <- df$genes$gene_id #this line assigns rownames to different tables in the df 
# STEP 5 look at the sample info in the df. A good idea is to compare with the sample info table
df$samples

# STEP 6 filter lowly expressed genes - run ?filterByExpr to see default arguments
keep <- filterByExpr(df) 
table(keep) #this line shows us how many genes are kept - TRUE
df <- df[keep, , keep.lib.sizes = FALSE] # this line actually filters the genes with low expression

# STEP 7 normalize the libraries and estimate dispersion (edgeR "equivalent" to computing means and sd for further statistics)
df <- calcNormFactors(df)
df <- estimateDisp(df)


# STEP 8 run the exact test suggested by edgeR to compare two groups
et <- exactTest(df, pair = c("ctrl", "rup")) # runs the test
et_results <- topTags(et, n = Inf) # save all (Inf) the results to et_results, instead of Inf you can put a number to only get top genes
et_results <- et_results$table %>% tibble() #convert the results table to data_frame (tibble is a newer tame for data_frame)

# STEP 9 see a summary of upregulated and downregulated genes
summary(decideTests(et))

# STEP 10 plot MDS - similar to a PCA plot (The multidimensional scaling (MDS) plot is frequently used to explore overall differences between samples)
plotMDS(df)

# STEP 11 plot MD (mean differences) - the genes in blue are downregulated, the genes in red upregulated
plotMD(et) 
abline(h = c(-1, 1), col = "blue") #adds to linest to the plot

# STEP 12 export the results to a file
et_results %>%
  write_tsv("/rnaseq-tutorials/data/with_ref/edger_results.csv")


