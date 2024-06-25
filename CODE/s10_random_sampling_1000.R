library(dplyr)

#### read in data ####
load("./results/top_sqtl_eqtl_group.RData")
load("./results/null_sqtl_group.RData") #5874411
load("./results/null_eqtl_group.RData") #11659741


#### Random Sampling ####
f_random <- function(dat_qtl, i, null_sqtl_group, null_eqtl_group){
  n_random <- nrow(dat_qtl[dat_qtl$group == i,])
  null_qtl <- c(null_sqtl_group[is.na(null_sqtl_group$group)==F & (null_sqtl_group$group == i), "SNP"], null_eqtl_group[(is.na(null_eqtl_group$group)==F) & (null_eqtl_group$group == i), "SNP"]) %>%
    unique()
  
  set.seed(2024)
  if (length(null_qtl)>=n_random){
    random_qtl <- data.frame(replicate(1000, sample(null_qtl, size=n_random, replace=FALSE)))
  } else {
    random_qtl <- data.frame(replicate(1000, sample(null_qtl, size=n_random, replace=TRUE)))
  }
}

## random sampling for top sQTLs
random_sqtl <- lapply(1:20, function(x) f_random(dat_sqtl, x, null_sqtl_group, null_eqtl_group))

## random sampling for top eQTLs
random_eqtl <- lapply(1:20, function(x) f_random(dat_eqtl, x, null_sqtl_group, null_eqtl_group))



#### save ####
save(random_sqtl, random_eqtl, file = "./results/sampling.RData")

