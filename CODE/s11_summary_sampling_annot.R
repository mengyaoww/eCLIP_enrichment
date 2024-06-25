#########################################################
#### SUMMARIZE ANNOTATION OF RANDOM SAMPLING RESULTS ####
#########################################################
library(dplyr)

#### read in data ####
load("./results/sampling.RData") ##1000 random sampling results
load("./results/annot_eqtlp4.RData") ##annotation of distinct cis-eQTL with p<10^(-4)
annot_eqtlp4 <- annot_nullsnp
load("./results/annot_nullsnp.RData") ##annotation of distinct cis-sQTL with p<10^(-4)
annot_sqtlp4 <- annot_nullsnp
#combine annotation of cis-eQTL and cis-sQTL
annot_nullsnp <- rbind(annot_eqtlp4, annot_sqtlp4) %>% distinct()



#### Summarize ramdom sampling for eQTL ####
random_eqtl_1 <- data.frame(unlist(lapply(random_eqtl, "[", , "X1")))
colnames(random_eqtl_1) <- "SNP"
annot_eqtl_1 <- left_join(random_eqtl_1, annot_nullsnp, by="SNP")
tab_eqtl_1 <- lapply(annot_eqtl_1[,4:211], table)
eqtl_rbp_sampling <- data.frame(lapply(tab_eqtl_1, function(x) data.frame(x)[1,2])) ## number of SNPs not annotated for each RBP

for (i in 2:1000){
  random_eqtl_rs <- data.frame(unlist(lapply(random_eqtl, "[", , paste0("X", i))))
  colnames(random_eqtl_rs) <- "SNP"
  annot_eqtl_rs <- left_join(random_eqtl_rs, annot_nullsnp, by="SNP")
  tab_eqtl_rs <- lapply(annot_eqtl_rs[,4:211], table)
  dat <- data.frame(lapply(tab_eqtl_rs, function(x) data.frame(x)[1,2]))
  
  eqtl_rbp_sampling <- rbind(eqtl_rbp_sampling, dat)
}

eqtl_rbp_sampling2 <- t(as.matrix(eqtl_rbp_sampling))
eqtl_rbp_sampling2 <- data.frame((15707 - eqtl_rbp_sampling2)/15707) ## proportion of SNPs annotated for each RBP
eqtl_rbp_sampling2$avg <- apply(eqtl_rbp_sampling2, 1, mean) ## average proportion of SNPs annotated for each RBP across 1000 samplings



#### Summarize random sampling for sQTL ####
random_sqtl_1 <- data.frame(unlist(lapply(random_sqtl, "[", , "X1")))
colnames(random_sqtl_1) <- "SNP"
annot_sqtl_1 <- left_join(random_sqtl_1, annot_nullsnp, by="SNP")
tab_sqtl_1 <- lapply(annot_sqtl_1[,4:211], table)
sqtl_rbp_sampling <- data.frame(lapply(tab_sqtl_1, function(x) data.frame(x)[1,2]))

for (i in 2:1000){
  random_sqtl_rs <- data.frame(unlist(lapply(random_sqtl, "[", , paste0("X", i))))
  colnames(random_sqtl_rs) <- "SNP"
  annot_sqtl_rs <- left_join(random_sqtl_rs, annot_nullsnp, by="SNP")
  tab_sqtl_rs <- lapply(annot_sqtl_rs[,4:211], table)
  dat2 <- data.frame(lapply(tab_sqtl_rs, function(x) data.frame(x)[1,2]))
  
  sqtl_rbp_sampling <- rbind(sqtl_rbp_sampling, dat2)
}

sqtl_rbp_sampling2 <- t(as.matrix(sqtl_rbp_sampling))
sqtl_rbp_sampling2 <- data.frame((4971 - sqtl_rbp_sampling2)/4971)
sqtl_rbp_sampling2$avg <- apply(sqtl_rbp_sampling2, 1, mean)



#### save ####
write.csv(eqtl_rbp_sampling2, file = "./results/eqtl_rbp_sampling.csv")
write.csv(sqtl_rbp_sampling2, file = "./results/sqtl_rbp_sampling.csv")





