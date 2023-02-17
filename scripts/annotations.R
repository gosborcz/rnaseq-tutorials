require(ggpubr)
require(tidyverse)
require(R.utils)
require(magrittr)


#this file does not contain many comments - we will read it together and figure out what the code means
#it is a good idea to add comments as we go through the file


setwd("/rnaseq-tutorials")

url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz"
dest <- "/rnaseq-tutorials/data/genecode.gtf.gz"

download.file(url,dest)
gunzip("/rnaseq-tutorials/data/genecode.gtf.gz")

gtf <- read_delim(
  "data/genecode.gtf")

# task 1 - look at the file. what is wrong with this table?
# add more parameters to read_delim above (don't worry about column X9 for now)

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

gtf[1:10,]

gtf <- gtf %>%
  filter(seqname == "chr12")

gtf[1:10,]

table(gtf$feature)

# task 2 - filter the gtf file to contain only transcripts



# insert the answer above this line


# now take a look at the column "attribute" - what is the problem?

att_names <- gtf %>%
  select(attribute) %>%
  apply(., MARGIN = 1, FUN = str_split, pattern = '"; ') %>%
  unlist() %>% trimws() %>% trimws(whitespace = ";") %>%
  sub(" .*$", "", .) %>% unique()

att_names
att_names <- att_names[att_names != ""]


for (att in att_names) {
  print(att)
}


for (att in att_names) {
  
  colatt <- apply(gtf, MARGIN = 1, function(x) {
    var <- str_extract(string = x[9],
                       pattern = sprintf('";\\s+%1$s[^;]+|^%1$s[^;]+;[^"]+"',att)) %>% 
      trimws(whitespace = '["; ]+', which = 'left') %>% 
      str_extract('(?<=")[^"]+(?=")')
    
  })
  gtf <- gtf %>% add_column("{att}" := colatt)
}

gtf %>% select(-c(attribute))

colnames(gtf)

# note: the code above is from here: https://www.biostars.org/p/272889/ and here: https://stackoverflow.com/questions/69901128/str-extract-regex-with-quotes-and-semicolons/69901344#69901344.
# Alternatively, there is a library rtracklayer with rtracklayer::import function that can import and parse a gtf file directly

transcripts_per_gene <- gtf$gene_name %>% 
  table() %>% 
  as.data.frame() %>% `colnames<-`(c("gene_name", "Freq"))

#task 3 - find 10 genes with the largest number of transcripts (hint: https://dplyr.tidyverse.org/reference/top_n.html)


gghistogram(transcripts_per_gene,
            x = "Freq",
            fill = "#00AFBB",
            add = "mean")

# task 4 - modify the histogram - transform the x axis into log2 (hint: http://www.sthda.com/english/wiki/ggplot2-axis-scales-and-transformations#axis-transformations)
# to add attributes to plots in ggplot2 a plus sign is used

#task 5 - make a table with transcript type frequency

#task 6* plot a bar plot with log10 axis, text rotated 45 degrees
# that will show frequencies of each transcript type



