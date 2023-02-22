# a trick to install the preprocessCore library in newer R verions:

BiocManager::install("preprocessCore", configure.args="--disable-threading")

# the script starts here 

require(plyr)
require(tidyverse)
require(preprocessCore)

# This script uses the de_novo data. Again, add your comments

files <- list.files(pattern = "quant.sf.genes", recursive = TRUE)

dfs <- ldply(
  files,
  function(dn) data.frame(Filename = dn, read_tsv(dn))
)

dfs$Filename <- sapply(
  strsplit(as.character(dfs$Filename), "/"), "[", 3
)


dfs <- pivot_wider(
  data = dfs,
  id_cols = c(Name),
  names_from = Filename,
  values_from = c("TPM")
)


annots <- read_tsv("data/de_novo/swissprot.blastx.outfmt6",
  col_names = FALSE
) %>%
  `colnames<-`(c(
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
  ))

annots$qseqid <- str_match(annots$qseqid, "(.*)_i(\\d)+$")[, 2]
annots$sseqid <- sapply(strsplit(as.character(annots$sseqid), "\\|"), "[", 2)
  
annots <- annots %>% distinct(
  qseqid,
  .keep_all = TRUE
)

# Task 1- The gene names we have are not "normal" gene names but uniprot names
# use the convert_to_gene_names.csv file in data/de_novo to annotate the annots file and then the dfs with human gene names


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


# your anwer here %>%
  p.adjust(method = 'fdr') %>%
  cbind(results, .) -> results

colnames(results)[10] <- "FDR"



# Task 6 calculate log2 fold-change (log2FC) for each detected transcript

# Task 7* - downolad data from this publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7691544/
# raw counts are here: http://149.156.177.112/projects/ifpan-kinga-dieta/cuffnorm/genes.count_table
# sample info here: http://149.156.177.112/projects/ifpan-kinga-dieta/cuffnorm/samples.table
# and here: http://149.156.177.112/projects/ifpan-kinga-dieta/analysis/samples-info.csv
#  and run your own two-factor analysis in a new file






