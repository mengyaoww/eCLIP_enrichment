###################################################
#### ANNOTATE NULL eQTLs with RBPs from ENCODE ####
###################################################
library(data.table)
library(dplyr)

n <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#### extract distinct cis-eQTL ####
#eqtl_p4 <- fread("https://ftp.ncbi.nlm.nih.gov/eqtl/original_submissions/phs002938_MolecularQTLs/eqtl/eqtl-result-combined-1e-4-annot-filtered.csv.gz")
#eqtl_p4_2 <- eqtl_p4 %>% filter(SNP_Chr == Tx_Chr) %>%
#  mutate(MAF = ifelse(is.na(EAF)==T, NA, ifelse(EAF <= 0.5, EAF, 1-EAF))) %>%
#  mutate(dis = abs(SNP_Pos - Tx_Start)) %>%
#  filter(dis <= 1e+6) ##238451941 --> 25822368 --> 11661223

#dat_eqtl <- eqtl_p4_2 %>% select(SNP, SNP_Chr, SNP_Pos) %>% distinct() #4335339
#save(dat_eqtl, file = "./results/distinct_eqtlp4.RData")


#### read in data ####
load("./results/distinct_eqtlp4.RData")
dat_eqtl$SNP_Chr <- paste0("chr", dat_eqtl$SNP_Chr)
if (n <= 433) {
  dat_snp <- dat_eqtl[((n-1)*10000+1):(n*10000),]
} else {
  dat_snp <- dat_eqtl[((n-1)*10000+1):nrow(dat_eqtl),]
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
    colnames(dat0)[ncol(dat0)] <- paste0(dateclip[dateclip$url==url, 2], "_", dateclip[dateclip$url==url, 4], "_", dateclip[dateclip$url==url, 5]) ##cellline_RBP_region
    
  }
  return(dat0)
}


system.time(null_snp_annot <- fmatch(dat_eclip, dat_snp))



#### save ####
save(null_snp_annot, file = paste0("./results/annot_eQTLp4/nulleqtl_annot_", n, ".RData"))

## each annotation file 10000 rows (SNPs) * 211 columns
## column 1: SNP
## column 2: SNP_Chr
## column 3: SNP_Pos
## column 4-211: annotation of each of 208 RBPs (Y: SNP is within RBP region; N: SNP is out of RBP region)


