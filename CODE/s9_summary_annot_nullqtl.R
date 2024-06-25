#################################################################
#### combine annotation outputs of null sQTLs and null eQTLs ####
#################################################################
library(dplyr)
library(purrr)

#### read in data ####
annot_sqtlp4 <- list.files("./annot_sQTLp4", full.names = TRUE, pattern = "\\.RData$") %>%
  map_df(~ get(load(file = .x)))
annot_eqtlp4 <- list.files("./annot_eQTLp4", full.names = TRUE, pattern = "\\.RData$") %>%
  map_df(~ get(load(file = .x)))


#### save ####
save(annot_sqtlp4, file = "./results/annot_sqtlp4.RData")
save(annot_eqtlp4, file = "./results/annot_eqtlp4.RData")
