###############################################################
#### ADDRESS EXPERIMANT INFORMATION downloaded from ENCODE ####
###############################################################
library(data.table)
library(dplyr)

setwd("/Users/mengyaow/Desktop/PROJECTS/sQTL/eCLIP")

##file contains eCLIP results from ENCODE; "*.bed.gz" files are default peak files (one for each experiment)
dat_eclip <- read.table("./ENCODE_file/eCLIP_url_208.txt", skip = 1)
eclip_list <- data.frame(dat_eclip[dat_eclip$V1 %like% "gz",])
colnames(eclip_list)[1] <- "url"
eclip_list$gzfile <- substr(eclip_list$url, 60,70)
##file contains information of 250 experiments from eCLIP ENCODE
dat_exp <- fread("./ENCODE_file/experiment_report_208.tsv", skip = 1, header = T)
dat_exp_e3 <- fread("./ENCODE_file/experiment_report_208_encode3.tsv", skip = 1, header = T)
dat_exp_e4 <- fread("./ENCODE_file/experiment_report_208_encode4.tsv", skip = 1, header = T)
dat_exp_e3$encode <- "3"
dat_exp_e3 <- dat_exp_e3[,c(1,41)]
dat_exp_e4$encode <- "4"
dat_exp_e4 <- dat_exp_e4[,c(1,41)]
dat_exp_e34 <- rbind(dat_exp_e3, dat_exp_e4)
dat_exp <- left_join(dat_exp, dat_exp_e34, by="ID")


## extract peak file from file column and remove unnecessary columns
a <- strsplit(dat_exp$Files, ",")

for (i in 1:nrow(dat_exp)){
  a2 <- lapply(a[[i]], function(x) substr(x,8, 18))
  dat_exp$peak_file[i] <- a2[lapply(a2, function(x) x%in%eclip_list$gzfile) == TRUE][[1]]
}
dat_exp2 <- dat_exp[,c("peak_file", "Biosample term name", "Target gene symbol", "Target of assay", "encode")]
dat_exp2 <- merge(dat_exp2, eclip_list, by.x = "peak_file", by.y = "gzfile")

## save
write.csv(dat_exp2, file = "./file/dat_experiment_208.csv", row.names = F)


