#############################################
#### SUMMARIZE ANNOTATED sQTLs and eQTLs ####
#############################################

#### eQTL ####
load("./results/match_eqtl_all_eclip208_sum.RData")
tab_eqtl <- lapply(eqtl_sum208[,21:228], table)
eqtl_freq <- data.frame(lapply(tab_eqtl, function(x) data.frame(x)[1,2]))
rbp <- colnames(eqtl_freq)
match <- as.numeric(eqtl_freq[1,]) #number of QTLs annotated with the corresponding RBP
eqtl_freq2 <- data.frame(eclip = rbp, n_match = match)
## 6 RBP have 15707 matches, they actually 0 match. We correct it.
eqtl_freq2[eqtl_freq2$n_match == 15707,]$n_match <- 0
write.csv(eqtl_freq2, file = "./results/annot_eqtl_all_freq.csv", row.names = F)



#### sQTL ####
load("./results/match_sqtl_all_eclip208_sum.RData")
tab_sqtl <- lapply(sqtl_sum208[,18:225], table)
sqtl_freq <- data.frame(lapply(tab_sqtl, function(x) data.frame(x)[1,2]))
rbp2 <- colnames(sqtl_freq)
match2 <- as.numeric(sqtl_freq[1,])
sqtl_freq2 <- data.frame(eclip = rbp2, n_match = match2)
## 12 RBP have 4971 matches, they actually 0 match. We correct it.
sqtl_freq2[sqtl_freq2$n_match == 4971,]$n_match <- 0
write.csv(sqtl_freq2, file = "./results/annot_sqtl_all_freq.csv", row.names = F)



