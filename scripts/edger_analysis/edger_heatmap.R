# STEP 0 Read libraries that are needed
library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(gplots)
require(preprocessCore)

# STEP 1 Read edger results (does not have counts in the table only results)
results <- read_tsv("/rnaseq-tutorials/data/with_ref/edger_results.csv")

# STEP 2 create vectors with group names and make a table with sample info
ctrls <- c("H27", "H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

sample_info <- cbind(
  c(ctrls,ruptured),
  c(rep("ctrl", 4),rep("rup",4)
  )
) %>% 
  as.data.frame() %>% `colnames<-`(c("samples", "groups")) 

# STEP 3 Read and transform the TPM counts for plotting
# the normalisation and file is the same as in the anova_analysis script
raw_data <- read_tsv(
  "/rnaseq-tutorials/data/with_ref/salmon.merged.gene_tpm.tsv"
)

tpms <- raw_data[,sample_info$samples]
tpms <- data.matrix(tpms)
tpms.norm <-  normalize.quantiles(tpms)
tpms.log <- log2(tpms.norm + 1)
colnames(tpms.log) <- sample_info$samples
rownames(tpms.log) <- raw_data$gene_id

# STEP 4 join the results table created above with the log-transformed counts
# we cannot simply use cbind as tmps.log have more rows than the results (these were filtered)
results <- cbind(
  results,
  tpms.log[match(results$gene_id, rownames(tpms.log)),]
  )


# STEP 5 apply filters for plotting - customise this as you wish e.g. add log2FC > 0 for upregulated genes only

to_plot <- results %>% 
  filter((FDR < 0.001) & (logCPM > 5)) %>% # this is an example filter for significantly regulated highly expressed genes
  select(sample_info$samples)

# STEP 6 add rownames to the to_plot table to show them on the heatmap
rownames(to_plot) <- results %>% 
  filter((FDR < 0.001) & (logCPM > 5)) %>%
  select(gene_name) %>% t()


# STEP 5 prepare color pallete for plotting
mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)


#STEP 7 prepare group names for plotting
group.names <- unique(as.character(sample_info$groups))
col.labels <- c(rep("", 2), group.names[1], rep(" ", 3),  group.names[2])

#STEP 8 create a function that cuts genes above certain threshold to a maximum value
cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}

#STEP 9 plot the heatmap - see heatmap.2 documentation for details https://www.rdocumentation.org/packages/gplots/versions/3.1.3/topics/heatmap.2
to_plot %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 2.5) %>%
  t %>%
  `colnames<-`(colnames(to_plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),
    trace="none",
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
