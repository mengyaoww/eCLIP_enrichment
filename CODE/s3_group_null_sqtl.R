##########################
#### group null sQTLs ####
##########################
library(dplyr)

#### read in data ####
load(paste0("./sqtl_batch12_selected", 1, ".Rdata"))
null_sqtl <- sqtl_batch12
for (i in 2:6){
  load(paste0("./sqtl_batch12_selected", i, ".Rdata"))
  null_sqtl <- rbind(null_sqtl, sqtl_batch12)
}



#### extract cis-sQTLs with distinct distance ####
null_sqtl2 <- null_sqtl %>% filter(SNP_Chr == Tx_Chr) %>%
  mutate(MAF = ifelse(is.na(EAF)==T, NA, ifelse(EAF <= 0.5, EAF, 1-EAF))) %>%
  mutate(dis = abs(SNP_Pos - Tx_Start)) %>%
  filter(dis <= 1e+6) #65620884 --> 10439098 --> 6227465
null_sqtl <- null_sqtl2 %>% select(SNP, SNP_Chr, SNP_Pos, MAF, dis) %>% distinct() #6227465 --> 5874411



#### group ####
null_sqtl_group <- null_sqtl %>%
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

table(null_sqtl_group$group, useNA = "always")



#### save ####
save(null_sqtl_group, file = "./results/null_sqtl_group.RData")


