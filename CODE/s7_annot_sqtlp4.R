###################################################
#### ANNOTATE NULL sQTLs with RBPs from ENCODE ####
###################################################
library(data.table)
library(dplyr)

n <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#### extract distinct cis-sQTLs ####
#load(paste0("./sqtl_batch12_selected", 1, ".Rdata"))
#null_sqtl <- sqtl_batch12
#for (i in 2:6){
#  load(paste0("./sqtl_batch12_selected", i, ".Rdata"))
#  null_sqtl <- rbind(null_sqtl, sqtl_batch12)
#}
#null_sqtl2 <- null_sqtl %>% filter(SNP_Chr == Tx_Chr) %>%
#  mutate(MAF = ifelse(is.na(EAF)==T, NA, ifelse(EAF <= 0.5, EAF, 1-EAF))) %>%
#  mutate(dis = abs(SNP_Pos - Tx_Start)) %>%
#  filter(dis <= 1e+6) #65620884 --> 10439098 --> 6227465
#dat_sqtl <- null_sqtl2 %>% select(SNP, SNP_Chr, SNP_Pos) %>% distinct() #6227465 --> 1814375
#save(dat_sqtl, file = "./results/distinct_sqtlp4.RData")


#### read in data ####
load("./results/distinct_sqtlp4.RData")
dat_sqtl$SNP_Chr <- paste0("chr", dat_sqtl$SNP_Chr)
if (n <= 181) {
  dat_snp <- dat_sqtl[((n-1)*10000+1):(n*10000),]
} else {
  dat_snp <- dat_sqtl[((n-1)*10000+1):nrow(dat_sqtl),]
}

dat_eclip <- read.csv("./file/dat_experiment_208.csv")



#### annotate ####
fmatch <- function(dateclip, dat0){
  for (j in 1:208){
    url <- dateclip$url[j]
    eclip <- fread(url, select = c(1:3))
    colnames(eclip) <- c("chr", "start", "end")
    eclip <- eclip[nchar(eclip$chr) <= 5,]
    
    for (i in 1:nrow(dat0)){
      a <- between(dat0$SNP_Pos[i], eclip[eclip$chr == dat0$SNP_Chr[i],]$start, eclip[eclip$chr == dat0$SNP_Chr[i],]$end)
      dat0$match[i] <- ifelse(all(a == FALSE), "N", "Y")
    }
    colnames(dat0)[ncol(dat0)] <- paste0(dateclip[dateclip$url==url, 2], "_", dateclip[dateclip$url==url, 4], "_", dateclip[dateclip$url==url, 5]) ##cellline_RBP_version
    
  }
  return(dat0)
}


system.time(null_snp_annot <- fmatch(dat_eclip, dat_snp))



#### save ####
save(null_snp_annot, file = paste0("./results/annot_sQTLp4/nullsnp_annot_", n, ".RData"))

## each annotation file 10000 rows (SNPs) * 211 columns
## column 1: SNP
## column 2: SNP_Chr
## column 3: SNP_Pos
## column 4-211: annotation of each of 208 RBPs (Y: SNP is within RBP region; N: SNP is out of RBP region)


