#######################################
#### group top eQTLs and top sQTLs ####
#######################################
library(dplyr)

#### read in data ####
dat_eqtl <- read.csv("./top_eqtl_plot.csv") #15707
dat_sqtl <- read.csv("./top_sqtl_plot.csv") #4971
dat_sqtl$s_MAF <- ifelse(dat_sqtl$s_EAF <= 0.5, dat_sqtl$s_EAF, 1-dat_sqtl$s_EAF)



#### group ####
match <- function(maf, tss){
  if (maf <= 0.05){
    if (tss <= 1e+3){
      a <- "1"
    } else if ((tss > 1e+3) & (tss <= 1e+4)){
      a <- "2"
    } else if ((tss > 1e+4) & (tss <= 1e+5)){
      a <- "3"
    } else {
      a <- "4"
    }
  } else if ((maf > 0.05) & (maf <= 0.1)){
    if (tss <= 1e+3){
      a <- "5"
    } else if ((tss > 1e+3) & (tss <= 1e+4)){
      a <- "6"
    } else if ((tss > 1e+4) & (tss <= 1e+5)){
      a <- "7"
    } else {
      a <- "8"
    }
  } else if ((maf > 0.1) & (maf <= 0.2)){
    if (tss <= 1e+3){
      a <- "9"
    } else if ((tss > 1e+3) & (tss <= 1e+4)){
      a <- "10"
    } else if ((tss > 1e+4) & (tss <= 1e+5)){
      a <- "11"
    } else {
      a <- "12"
    } 
  } else if ((maf > 0.2) & (maf <= 0.3)){
    if (tss <= 1e+3){
      a <- "13"
    } else if ((tss > 1e+3) & (tss <= 1e+4)){
      a <- "14"
    } else if ((tss > 1e+4) & (tss <= 1e+5)){
      a <- "15"
    } else {
      a <- "16"
    }
  } else {
    if (tss <= 1e+3){
      a <- "17"
    } else if ((tss > 1e+3) & (tss <= 1e+4)){
      a <- "18"
    } else if ((tss > 1e+4) & (tss <= 1e+5)){
      a <- "19"
    } else {
      a <- "20"
    }
  }
  return(a)
}

group_eqtl <- lapply(1:15707, function(x) match(dat_eqtl$MAF[x], dat_eqtl$dis[x]))
for (i in 1:15707){
  dat_eqtl$group[i] <- group_eqtl[[i]]
}
table(dat_eqtl$group, useNA = "always")

group_sqtl <- lapply(1:4971, function(x) match(dat_sqtl$s_MAF[x], dat_sqtl$dis[x]))
for (i in 1:4971){
  dat_sqtl$group[i] <- group_sqtl[[i]]
}
table(dat_sqtl$group, useNA = "always")


save(dat_sqtl, dat_eqtl, file = "./results/top_sqtl_eqtl_group.RData")

