##############################
#### 1000 RANDOM SAMPLING ####
##############################
library(dplyr)
library(ggplot2)

#### read ind data ####
dat_eclip <- read.csv("./file/dat_experiment_208.csv")
## find duplicated rbp in eCLIP results 
## (There are 5 RBPs have duplicated experiments with same cell line between version3 and version4, I use results from version4)
ind <- dat_eclip[,c(2,4)] 
dat_dup <- dat_eclip[duplicated(ind) | duplicated(ind, fromLast=TRUE),] %>%
  mutate(rbp = paste0(Biosample.term.name, "_", Target.of.assay, "_", encode)) %>% filter(encode == 3)

annot_sqtl <- read.csv("./annot_sqtl_all_freq.csv") %>%
  filter(!eclip %in% dat_dup$rbp)
colnames(annot_sqtl)[2] <- "match_sqtl"
annot_sqtl$sqtl_p <- (annot_sqtl$match_sqtl)/4971
annot_eqtl <- read.csv("./results/annot_eqtl_all_freq.csv") %>%
  filter(!eclip %in% dat_dup$rbp)
colnames(annot_eqtl)[2] <- "match_eqtl"
annot_eqtl$eqtl_p <- (annot_eqtl$match_eqtl)/15707

sqtl_rbp_sampling <- read.csv("./results/sqtl_rbp_sampling.csv") %>%
  filter(!X %in% dat_dup$rbp)
colnames(sqtl_rbp_sampling)[1] <- "eclip"
enrich_sqtl <- left_join(annot_sqtl, sqtl_rbp_sampling[,c(1,1002)], by = "eclip") %>% mutate(fold = sqtl_p/avg)

eqtl_rbp_sampling <- read.csv("./results/eqtl_rbp_sampling.csv") %>%
  filter(!X %in% dat_dup$rbp)
colnames(eqtl_rbp_sampling)[1] <- "eclip"
enrich_eqtl <- left_join(annot_eqtl, eqtl_rbp_sampling[,c(1,1002)], by = "eclip") %>% mutate(fold = eqtl_p/avg)



#### violin plot ####
identical(enrich_eqtl$eclip, enrich_sqtl$eclip) #TRUE
annot_all <- data.frame(
  eclip = rep(enrich_eqtl$eclip, 2),
  enrich = c(enrich_sqtl$fold, enrich_eqtl$fold),
  qtl = rep(c("sQTL", "eQTL"), each = 203)
)
mean_annot <- data.frame(
  qtl = c("sQTL", "eQTL"),
  mean = paste0(c(formatC(mean(enrich_sqtl$fold), digits = 1, format = "f"), 
                  formatC(mean(enrich_eqtl$fold), digits = 1, format = "f")), "-fold")
)


v1 <- ggplot(annot_all, aes(x = qtl, y = enrich)) +
  geom_violin(aes(fill = qtl)) +
  geom_boxplot(width=0.03, fill = "white") +
  geom_text(data = mean_annot, aes(x = qtl, y = 6, label = mean, color = qtl), hjust = -0.5) +
  labs(x = "", y = "Fold enrichment", title = "eCLIP peaks of RBPs (ENCODE)") +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title.y=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none")

v2 <- ggplot(annot_all, aes(x = qtl, y = enrich)) +
  geom_violin(aes(fill = qtl)) +
  geom_boxplot(width=0.03, fill = "white") +
  geom_text(data = mean_annot, aes(x = qtl, y = 6, label = mean, color = qtl), hjust = -0.5) +
  labs(x = "", y = "Fold enrichment", title = "eCLIP peaks of RBPs (ENCODE)") +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title.y=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1))

v3 <- ggplot(annot_all, aes(x = qtl, y = enrich)) +
  geom_violin(aes(fill = qtl)) +
  geom_boxplot(width=0.03, fill = "white") +
  geom_text(data = mean_annot, aes(x = qtl, y = 6, label = mean, color = qtl), hjust = -0.5) +
  labs(x = "", y = "Fold enrichment", title = "eCLIP peaks of RBPs (ENCODE)") +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title.y=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none")+
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

v4 <- ggplot(annot_all, aes(x = qtl, y = enrich)) +
  geom_violin(aes(fill = qtl)) +
  geom_boxplot(width=0.03, fill = "white") +
  geom_text(data = mean_annot, aes(x = qtl, y = 6, label = mean, color = qtl), hjust = -0.5) +
  labs(x = "", y = "Fold enrichment", title = "eCLIP peaks of RBPs (ENCODE)") +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title.y=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))



#### save ####
pdf(file = "./results/plot_eclip.pdf", width = 6, height = 4)
v1
v2
v3
v4
dev.off()



