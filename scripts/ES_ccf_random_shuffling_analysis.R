# Load unfiltered mutation data
all_mut_ccf <- readRDS("/home/labgroups/ccgg/downloads/pub/wu_2022/Immunoediting/data/all_mut_ccf_tpm.rds")

# Remove data without CCF information
all_mut_ccf<- all_mut_ccf[!is.na(all_mut_ccf$ccf_hat),]

# Shuffle neoantigen labelling & CCF
sim_ccf_df<- data.frame(
  sample=all_mut_ccf$sample,
  neo=sample(all_mut_ccf$neo),
  ccf=all_mut_ccf$ccf_hat
)

# Filter at least 1 neoantigenic and 1 subclonal mutation (CCF<0.6)
neo_t<- table(sim_ccf_df$sample, sim_ccf_df$neo)
samples_incl_neo<- rownames(neo_t[neo_t[,"neo"]>0,]) 
ccf_t<- table(sim_ccf_df$sample, sim_ccf_df$ccf<0.6)
samples_incl_clon<- rownames(ccf_t[ccf_t[,"TRUE"]>0,]) 
sim_ccf_df_filt<- sim_ccf_df[sim_ccf_df$sample%in%intersect(samples_incl_neo,samples_incl_clon),]

dt<- data.frame(
  sample=sim_ccf_df_filt$sample,
  neo="no",
  ccf=sim_ccf_df_filt$ccf
)
dt$neo[sim_ccf_df_filt$neo=="neo"]<- "yes"

# # Save
# saveRDS(dt,file = "temp/mut_df_ccf.rds")

#################################################
# Calculate ES RNA & plot
#################################################

dt<- readRDS(file = "temp/mut_df_ccf.rds")

samples<- unique(dt$sample)
source("scripts/functions/cal_es_new_test.R")
es_ccf<- mclapply(samples, function(x) cal_es_new_test(dt[dt$sample==x,]),mc.cores = 10)
summary(as.numeric(es_ccf)) # -01434 
names(es_ccf)<- substr(samples,1,15)

# Median ES per cancer
load("../immunoediting_2019/data/TCGA_manifest.RData")
ES_med<- sort(tapply(as.numeric(es_ccf), TCGA_cancer_id[names(es_ccf),1], "median", na.rm=T))

# Calculate p values from permutation data
# see scripts/ES_ccf_perm.R

# Calculate p values from permutation data
perm_matrix<- readRDS(file = "temp/perm_matrix_ccf.rds")
perm_matrix_cancer<- t(apply(perm_matrix, 1, function(x) tapply(x, TCGA_cancer_id[substr(colnames(perm_matrix),1,15),1], "median",na.rm=T)))

pval_can<- rep(NA, ncol(perm_matrix_cancer))
names(pval_can)<- colnames(perm_matrix_cancer)
for(c in colnames(perm_matrix_cancer)){
  pval_can[c]<- mean(ES_med[c]>perm_matrix_cancer[,c])
}
library(matrixStats)
pval_pan<- mean(median(as.numeric(es_ccf),na.rm=T)>rowMedians(perm_matrix, na.rm=T))

# plot
es_df<- data.frame(es=as.numeric(es_ccf), cancer=TCGA_cancer_id[names(es_ccf),1], sample=names(es_ccf))
es_df$cancer<- factor(es_df$cancer, levels=names(ES_med))
library(ggprism)
psign_can<- rep("ns",length(pval_can))
psign_can[pval_can<0.05]<- "*"
psign_can[pval_can<0.01]<- "**"
psign_can[pval_can<0.001]<- "***"
pval_can_df<- data.frame(p=pval_can, sign=psign_can)
pval_can_df<- pval_can_df[levels(es_df$cancer),]
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
    x=paste0("Median = ",round(median(es_df$es, na.rm=T),digits = 3),"\n n = ",nrow(es_df)),y="ES CCF",size=1)+
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

# demonstrate remaining TMB correlations
#########################################

# By sample
all_mut_ccf <- readRDS("/home/labgroups/ccgg/downloads/pub/wu_2022/Immunoediting/data/all_mut_ccf_tpm.rds")
TMB<- table(all_mut_ccf$sample) # Original, before filtering

neo_t<- table(dt$sample, dt$neo)
prop_neo<- prop.table(neo_t,1)[,"yes"]

sim_ccf_df_sample<- data.frame(
  sample=samples,
  es=unlist(es_ccf),
  tmb=as.numeric(TMB[samples]),
  neo=prop_neo[samples]
)
load("../immunoediting_2019/data/TCGA_manifest.RData")
sim_ccf_df_sample$cancer<- TCGA_cancer_id[substr(sim_ccf_df_sample$sample,1,15),1] 

# By cancer
es_cancer<- sort(tapply(sim_ccf_df_sample$es, sim_ccf_df_sample$cancer, "median", na.rm=T))
tmb_cancer<- sort(tapply(sim_ccf_df_sample$tmb, sim_ccf_df_sample$cancer, "median", na.rm=T))
neo_cancer<- sort(tapply(sim_ccf_df_sample$neo, sim_ccf_df_sample$cancer, "median", na.rm=T))

cancers<- names(tmb_cancer)
cancer_df<- data.frame(
  cancer=cancers,
  es=es_cancer[cancers],
  tmb=tmb_cancer[cancers],
  neo=neo_cancer[cancers]
)

library(ggrepel)
library(ggpubr)
p_tmb_es<- ggplot(cancer_df, aes(x=tmb, y=es, label=cancer)) +
  geom_point() +
  geom_text_repel(size=2) +
  theme_prism()+
  xlab("TMB") +
  ylab("ES CCF") +
  theme(
    axis.line = element_line(size=0.4),
    axis.ticks = element_line(size=0.4),
    text = element_text(size = 7)
  ) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  stat_cor(method = "spearman", label.x = 200, hjust=1, size=2)

p_tmb_neo<- ggplot(cancer_df, aes(x=tmb, y=neo, label=cancer)) +
  geom_point() +
  geom_text_repel(size=2) +
  scale_y_continuous(labels = scales::percent) +
  theme_prism()+
  xlab("TMB") +
  ylab("Neoantigenic mutations (%)") +
  theme(
    axis.line = element_line(size=0.4),
    axis.ticks = element_line(size=0.4),
    text = element_text(size = 7)
  ) +
  stat_cor(method = "spearman", label.x = 200, hjust=1, size=2)

# Filter effect?
samples<- unique(all_mut_ccf$sample)
n_samples_preFilter<- table(TCGA_cancer_id[substr(samples,1,15),1])
n_samples_postFilter<- table(sim_ccf_df_sample$cancer)
pct_postfilter<- n_samples_postFilter[names(n_samples_preFilter)]/n_samples_preFilter

cancer_df$pct_filter<- pct_postfilter[cancer_df$cancer]

p_tmb_filter<- ggplot(cancer_df, aes(x=tmb, y=pct_filter, label=cancer)) +
  geom_point() +
  geom_text_repel(size=2) +
  scale_y_continuous(labels = scales::percent) +
  theme_prism()+
  xlab("TMB") +
  ylab("Included samples after filtering (%)") +
  theme(
    axis.line = element_line(size=0.4),
    axis.ticks = element_line(size=0.4),
    text = element_text(size = 7)
  ) +
  stat_cor(method = "spearman", label.x = 200, label.y = floor(10*min(cancer_df$pct_filter))/10, hjust=1, vjust=0, size=2)

# Merge plots
p_cancer_tmb<- plot_grid(
  p_tmb_es, p_tmb_neo, p_tmb_filter,
  ncol=3,
  labels=c("F","G","H")
)

p<- plot_grid(
  p, p_cancer_tmb,
  ncol=1,
  labels=c("E",NA)
)

# Save
########
saveRDS(p, "results/data/p_shuffle_ccf.rds")





