### this is a tutorial on how to parse a gtf file ###
### remember to change the file paths to your own

# STEP 1 load required libraries
require(ggpubr)
require(tidyverse)
require(R.utils)
require(magrittr)

# STEP 2 set working directory
setwd("/rnaseq-tutorials")

# STEP 3 save url of gtf as src_file variable and destination path as dest
src_file <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz"
dest <- "/rnaseq-tutorials/data/genecode.gtf.gz"

# STEP 4 download gtf file
download.file(src_file,dest)

# STEP 5 unzip gft file
gunzip("/rnaseq-tutorials/data/genecode.gtf.gz")

# STEP 6 read the gtf file in ignoring the lines starting with ##
gtf <- read_delim(
  "data/genecode.gtf",
  delim="\t",
  comment = "##",
  col_names = FALSE)

# STEP 7 rename columns ing gtf
colnames(gtf) <- c("seqname",
                   "source",
                   "feature",
                   "start",
                   "end",
                   "score",
                   "strand",
                   "frame",
                   "attribute"
                   )
# STEP 8 filter gtf to chr12 only
gtf <- gtf %>%
  filter(seqname == "chr12")

# STEP 9 show numbers of each different feature
table(gtf$feature)

# STEP 10 fiter gtf to contain only transcript features
gtf <- gtf %>% 
  filter(feature == "transcript")


# STEP 11 the script below is from biostars forum  https://www.biostars.org/p/272889/ 
# it parses the attribute column into separate columns with different names

att_names <- gtf %>%
  select(attribute) %>%
  apply(., MARGIN = 1, FUN = str_split, pattern = '"; ') %>%
  unlist() %>% trimws() %>% trimws(whitespace = ";") %>%
  sub(" .*$", "", .) %>% unique()

att_names <- att_names[att_names != ""]


for (att in att_names) {
  
  colatt <- apply(gtf, MARGIN = 1, function(x) {
    var <- str_extract(string = x[9],
                       pattern = sprintf('";\\s+%1$s[^;]+|^%1$s[^;]+;[^"]+"',att)) %>% 
      trimws(whitespace = '["; ]+', which = 'left') %>% 
      str_extract('(?<=")[^"]+(?=")')
    
  })
  gtf <- gtf %>% add_column("{att}" := colatt)
}

gtf <- gtf %>% 
  select(-attribute)

# STEP 12 look at the colnames of the gtf filre after running the script
colnames(gtf)

# STEP 13 make a new data frame with transcripts per gene
transcripts_per_gene <- gtf$gene_name %>% 
  table() %>% 
  as.data.frame() %>% `colnames<-`(c("gene_name", "count"))

# STEP 14 show 10 genes with the highest number of transcripts
transcripts_per_gene %>%
  top_n(10) #number can be replaced with the requested number

# STEP 15 make a histogram of the number of genes (see documentation of ggpubr for more details)
gghistogram(transcripts_per_gene,
            x = "count",
            fill = "#00AFBB",
            add = "mean") + scale_x_continuous(trans='log2')

# STEP 16 make a table with the numbers of each transcript type
transcript_types <- gtf$transcript_type %>% 
  table() %>% 
  as.data.frame() %>%
  `colnames<-`(c("transcript_type", "transcript_count")) %>%
  arrange(transcript_count) 

# STEP 16 make a barplot with different transcript types (see documentation of ggpubr for more details)
ggbarplot(transcript_types, "transcript_type", "transcript_count",
          fill = "steelblue", color = "steelblue",
          order = transcript_types$transcript_type,
          label = TRUE, lab.pos = "in", lab.col = "white") + scale_y_continuous(
            trans='log10'
          ) + rotate_x_text(angle = 45)



