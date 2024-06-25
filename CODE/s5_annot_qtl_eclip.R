######################################################
#### MATCH eQTL RESULTS WITH ENCODE eCLIP RESULTS ####
######################################################
library(data.table)

#### read in data ####
dat_sqtl_all <- read.csv("top_sqtl_plot.csv") #4971 sQTLs
dat_sqtl_all$s_SNP_Chr <- paste0("chr", dat_sqtl_all$s_SNP_Chr)
dat_eqtl_all <- read.csv("top_eqtl_plot.csv") #15707 eQTLs
dat_eqtl_all$SNP_Chr <- paste0("chr", dat_eqtl_all$SNP_Chr)
## addressed file contains eCLIP results from ENCODE; "*.bed.gz" files are default peak files (one for each experiment)
dat_eclip <- read.csv("./file/dat_experiment_208.csv")



#### match eQTL position and eCLIP regions ####
match_eqtl_eclip <- function(qtl, url){
  eclip <- fread(url)
  eclip <- eclip[,1:8]
  colnames(eclip) <- c("chr", "start", "end", "label", "1000", "strand", "log2(fold)", "-log10(P)")
  eclip <- eclip[nchar(eclip$chr) <= 5,] #limit the chromosome to chr1-22,X,Y,M
  for (i in 1:nrow(qtl)){
    a <- between(qtl$SNP_Pos[i], eclip[eclip$chr == qtl$SNP_Chr[i],]$start, eclip[eclip$chr == qtl$SNP_Chr[i],]$end)
    qtl$match[i] <- ifelse(all(a == FALSE), "N", paste0(dat_eclip[dat_eclip$url==url, 2], "_", dat_eclip[dat_eclip$url==url, 4]))
    qtl$n_match[i] <- length(a[a == TRUE])
  }
  colnames(qtl)[21] <- paste0(dat_eclip[dat_eclip$url==url, 2], "_", dat_eclip[dat_eclip$url==url, 4], "_", dat_eclip[dat_eclip$url==url, 5])
  return(qtl)
}

eqtl_list208 <- lapply(dat_eclip$url[1:208], function(x) match_eqtl_eclip(dat_eqtl_all, x))



#### match sQTL position and eCLIP regions ####
match_sqtl_eclip <- function(qtl, url){
  eclip <- fread(url)
  eclip <- eclip[,1:8]
  colnames(eclip) <- c("chr", "start", "end", "label", "1000", "strand", "log2(fold)", "-log10(P)")
  eclip <- eclip[nchar(eclip$chr) <= 5,] #limit the chromosome to chr1-22,X,Y,M
  for (i in 1:nrow(qtl)){
    a <- between(qtl$s_SNP_Pos[i], eclip[eclip$chr == qtl$s_SNP_Chr[i],]$start, eclip[eclip$chr == qtl$s_SNP_Chr[i],]$end)
    qtl$match[i] <- ifelse(all(a == FALSE), "N", paste0(dat_eclip[dat_eclip$url==url, 2], "_", dat_eclip[dat_eclip$url==url, 4]))
    qtl$n_match[i] <- length(a[a == TRUE])
  }
  colnames(qtl)[18] <- paste0(dat_eclip[dat_eclip$url==url, 2], "_", dat_eclip[dat_eclip$url==url, 4], "_", dat_eclip[dat_eclip$url==url, 5])
  return(qtl)
}

sqtl_list208 <- lapply(dat_eclip$url[1:208], function(x) match_sqtl_eclip(dat_sqtl_all, x))



#### summarize matching results ####
library(dplyr)
eqtl_sum208 <- eqtl_list208[[1]]
eqtl_sum208 <- eqtl_sum208[-22]
for (i in 2:208){
  eqtl_sum208 <- bind_cols(eqtl_sum208, eqtl_list208[[i]][21])
}

sqtl_sum208 <- sqtl_list208[[1]]
sqtl_sum208 <- sqtl_sum208[-19]
for (i in 2:208){
  sqtl_sum208 <- bind_cols(sqtl_sum208, sqtl_list208[[i]][18])
}

#first 20 columns: information about each SNP
#column 21- 228: annotation of SNPs for each RBP (Y: SNP is within the corresponding RBP regions; N: SNP is not within the corresponding RBP regions)



#### save ####
save(eqtl_sum208, file = "./results/match_eqtl_all_eclip208_sum.RData")
save(sqtl_sum208, file = "./results/match_sqtl_all_eclip208_sum.RData")



##if we use all peak files from ENCODE, we removed experiments of K562/HepG2 with same target 
##(duplicated experiments between ENCODE3 and ENCODE4)

#### WE REMOVE REGIONS WITH "chrUn" OR chrN_random" OR "chrUn_random"
## ChrUn contains clone contigs that cannot be confidently placed on a specific chromosome. 
## For the chrN_random and chrUn_random files, we essentially just concatenate together all the contigs into short pseudo-chromosomes. 
## The coordinates of these are fairly arbitrary, although the relative positions of the coordinates are good within a contig


