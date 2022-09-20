# 1) Somatic mutation data Wu et al.
######################################

all_mut_exp <- readRDS("downloads/pub/wu_2022/Immunoediting/data/all_mut_tpm_not_filter.rds")

# 2) Calculate ES
###################

source("scripts/functions/cales.R")
source("scripts/functions/cales_t.R")

# Create mutation datatable as specified by readme in https://github.com/wt12318/NeoEnrichment/
mut_dt<- data.frame(
  gene = all_mut_exp$gene,
  sample = substr(all_mut_exp$sample,1,15),
  neo = all_mut_exp$neo,
  exp = all_mut_exp$tpm_exp
)

# Calculate ES RNA
library(tidyverse)
samples<- unique(mut_dt$sample)
ES_rna<- lapply(samples, function(x) cales_t(data = mut_dt, barcode = x, calp = F, cal_type = "exp", type = "II", sample_counts = 1000))
summary(as.numeric(ES_rna)) # -0.041
names(ES_rna)<- samples

# 3. Permutation test
#####################
saveRDS(mut_dt,file = "temp/mut_df_wu.rds")
# see scripts/ES_rna_wu_perm.R

# 4. Plot
#########

# Median ES per cancer
load("../immunoediting_2019/data/TCGA_manifest.RData")
ES_med<- sort(tapply(as.numeric(ES_rna), TCGA_cancer_id[names(ES_rna),1], "median", na.rm=T))

# Calculate p values from permutation data
perm_matrix<- readRDS(file = "temp/perm_matrix_wu.rds")
perm_matrix_cancer<- t(apply(perm_matrix, 1, function(x) tapply(x, TCGA_cancer_id[colnames(perm_matrix),1], "median",na.rm=T)))

pval_can<- rep(NA, ncol(perm_matrix_cancer))
names(pval_can)<- colnames(perm_matrix_cancer)
for(c in colnames(perm_matrix_cancer)){
  pval_can[c]<- mean(ES_med[c]>perm_matrix_cancer[,c])
}
pval_pan<- mean(median(as.numeric(ES_rna),na.rm=T)>rowMedians(perm_matrix, na.rm=T))

# plot
es_df<- data.frame(es=as.numeric(ES_rna), cancer=TCGA_cancer_id[names(ES_rna),1], sample=names(ES_rna))
es_df$cancer<- factor(es_df$cancer, levels=names(ES_med))
library(ggprism)
psign_can<- rep("ns",length(pval_can))
psign_can[pval_can<0.05]<- "*"
psign_can[pval_can<0.01]<- "**"
psign_can[pval_can<0.001]<- "***"
pval_can_df<- data.frame(p=pval_can, sign=psign_can)
pval_can_df$x<- 1:nrow(pval_can_df)
pval_can_df$y<- 1.1
pval_can_df$size<- 3
pval_can_df$size[pval_can_df$sign=="ns"]<- 2
pval_can_df$adj<- 1
pval_can_df$adj[pval_can_df$sign=="ns"]<- 0.7

p_pan<- ggplot(data=es_df,aes(x=1,y=es))+
  geom_violin(alpha=0.7,width=0.5)+
  geom_boxplot(width=0.2, outlier.size = 0.5)+
  theme_prism()+
  labs(
    x=paste0("Median = ",round(median(es_df$es, na.rm=T),digits = 3),"\n n = ",nrow(es_df)),y="ES RNA",size=1)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(size=0.4),
        axis.ticks = element_line(size=0.4)
  )+
  theme(text = element_text(size = 7)) +
  annotate(geom="text", x=1, y=1.1, label=paste0("p = ",signif(pval_pan, digits=1)),
           color="red",size=3)

p_can<- ggplot(data=es_df,aes(x=cancer,y=es))+
  geom_boxplot(outlier.size = 0.5)+
  theme_prism()+
  theme(axis.title.x = element_blank(),
        axis.line = element_line(size=0.4),
        axis.ticks = element_line(size=0.4)
  )+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1))+
  theme(text = element_text(size = 6))+
  ylab("")+
  geom_text(data=pval_can_df,aes(x=x,y=y,label = sign),size=pval_can_df$size, vjust=pval_can_df$adj)+
  geom_hline(yintercept=0,
             color = "red", size=1)

p<- plot_grid(
  p_pan, p_can,
  ncol = 2,
  rel_widths = c(1,4),
  axis = "rlbt",
  align = "hv")

# save plot
saveRDS(p, file = "results/data/p_wu.rds")
ggsave("results/figs/fig1a.pdf", p, width = 178, height = 265/3, units = "mm")

# Save results
##############
save.image(file = "results/data/results_wu.RData")
