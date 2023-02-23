require(magrittr)

setwd("/seq/intelliseq/malus_domestica")
sample_info <- read.csv("/seq/intelliseq/malus_domestica/SraRunInfo.csv", header = TRUE, colClasses = "character")

get_sample_info <- function(x) {
  split <- unlist(strsplit(x, "_"))
  variety = split[1]
  if (length(split) > 2) {
    treat = "ctrl"
    time = split[3]
  } else {
    treat = "blight"
    time = split[2]
  }
  out <- c(variety, treat, time)
  names(out) <- c("variety", "treat", "time")
  out
}

sample_info$SampleName %>% lapply(get_sample_info) %>% do.call(rbind, .) %>% data.frame
