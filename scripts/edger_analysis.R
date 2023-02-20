#installing EdgeR package - run this only once

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

#script starts here, remember to add your comments

require(edgeR)
require(tidyverse)

#Task 1 - read the /rnaseq-tutorials/data/with_ref/salmon.merged.gene_counts.tsv file into R

rawdata <- #put your answer here
 
y <- DGEList(counts=rawdata[,3:10],
             genes=rawdata[,1:2])

rownames(y$counts) <- rownames(y$genes) <- y$genes$gene_name

keep <- filterByExpr(y)

table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

plotMDS(y)

ctrls <- c("H27","H45", "H91", "H74")
ruptured <- c("H2", "H47", "H56", "H87")

y$samples$group <- as.factor(c("rup",
                               "ctrl",
                               "ctrl",
                               "rup",
                               "rup",
                               "ctrl",
                               "rup",
                               "ctrl"))

y$samples

y <- estimateDisp(y)
et <- exactTest(y, pair=c("ctrl","rup"))
et_results <- topTags(et, n=Inf)
et_results <- et_results$table %>% tibble()

summary(decideTests(et))
plotMD(et)
abline(h=c(-1, 1), col="blue")

#Task 2 - annotate the results with transcript types from the gtf file

# Task 3 - use enrichR (manually or via R) to look for overrepresented pathways in the results

# Task 4* - make a new script edger_analysis_de_novo.R where you will analyse
# the de novo assembly with edger [take a look at the anova_analysis for hints]

# Task 5* - downolad data from this publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7691544/
# raw counts are here: http://149.156.177.112/projects/ifpan-kinga-dieta/cuffnorm/genes.count_table
# sample info here: http://149.156.177.112/projects/ifpan-kinga-dieta/cuffnorm/samples.table
# and here: http://149.156.177.112/projects/ifpan-kinga-dieta/analysis/samples-info.csv
#  and run your own two-factor analysis in a new file

