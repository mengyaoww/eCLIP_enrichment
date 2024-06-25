##########################
#### group null eQTLs ####
##########################
library(dplyr)
library(data.table)

#### read in data ####
eqtl_p4 <- fread("https://ftp.ncbi.nlm.nih.gov/eqtl/original_submissions/phs002938_MolecularQTLs/eqtl/eqtl-result-combined-1e-4-annot-filtered.csv.gz")



#### extract cis-eQTLs with distinct distance ####
eqtl_p4_2 <- eqtl_p4 %>% filter(SNP_Chr == Tx_Chr) %>%
  mutate(MAF = ifelse(is.na(EAF)==T, NA, ifelse(EAF <= 0.5, EAF, 1-EAF))) %>%
  mutate(dis = abs(SNP_Pos - Tx_Start)) %>%
  filter(dis <= 1e+6) ##238451941 --> 25822368 --> 11661223
null_eqtl <- eqtl_p4_2 %>% select(SNP, SNP_Chr, SNP_Pos, MAF, dis) %>% distinct() #11659741



#### group ####
null_eqtl_group <- null_eqtl %>%
  mutate(group = ifelse(MAF <= 0.05 & dis <= 1e+3, "1", 
                        ifelse(MAF <= 0.05 & dis>1e+3 & dis <= 1e+4, "2",
                               ifelse(MAF <= 0.05 & dis>1e+4 & dis <= 1e+5, "3",
                                      ifelse(MAF <= 0.05 & dis>1e+5 & dis <= 1e+6, "4", 
                                             ifelse(MAF > 0.05 & MAF <= 0.1 & dis <= 1e+3, "5",
                                                    ifelse(MAF > 0.05 & MAF <= 0.1 & dis>1e+3 & dis <= 1e+4, "6",
                                                           ifelse(MAF > 0.05 & MAF <= 0.1 & dis>1e+4 & dis <= 1e+5, "7",
                                                                  ifelse(MAF > 0.05 & MAF <= 0.1 & dis>1e+5 & dis <= 1e+6, "8",
                                                                         ifelse(MAF > 0.1 & MAF <= 0.2 & dis <= 1e+3, "9",
                                                                                ifelse(MAF > 0.1 & MAF <= 0.2 & dis>1e+3 & dis <= 1e+4, "10",
                                                                                       ifelse(MAF > 0.1 & MAF <= 0.2 & dis>1e+4 & dis <= 1e+5, "11",
                                                                                              ifelse(MAF > 0.1 & MAF <= 0.2 & dis>1e+5 & dis <= 1e+6, "12",
                                                                                                     ifelse(MAF > 0.2 & MAF <= 0.3 & dis <= 1e+3, "13",
                                                                                                            ifelse(MAF > 0.2 & MAF <= 0.3 & dis>1e+3 & dis <= 1e+4, "14",
                                                                                                                   ifelse(MAF > 0.2 & MAF <= 0.3 & dis>1e+4 & dis <= 1e+5, "15",
                                                                                                                          ifelse(MAF > 0.2 & MAF <= 0.3 & dis>1e+5 & dis <= 1e+6, "16",
                                                                                                                                 ifelse(MAF > 0.3 & dis <= 1e+3, "17",
                                                                                                                                        ifelse(MAF > 0.3 & dis>1e+3 & dis <= 1e+4, "18",
                                                                                                                                               ifelse(MAF > 0.3 & dis>1e+4 & dis <= 1e+5, "19", "20"))))))))))))))))))))

table(null_eqtl_group$group, useNA = "always")

check1 <- null_eqtl[is.na(null_eqtl$MAF)==T,] #MAF NA here are also EAF NA in original eQTL results datatset



#### save ####
save(null_eqtl_group, file = "./results/null_eqtl_group.RData")



