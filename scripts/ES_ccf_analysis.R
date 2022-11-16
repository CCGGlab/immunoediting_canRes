##################
# ES CCF analysis
##################

library(ggplot2)
library(ggprism)
library(cowplot)
library(dplyr)
library(ggrepel)
library(ggpubr)

# Check correlation to TMB
############################

# Load ES CCF data used in figure 2A
es_ccf <- readRDS("/home/labgroups/ccgg/downloads/pub/wu_2022/Immunoediting/data/neo_nes_ccf06_1.rds")

# Get TMB from mutation data
all_mut_ccf <- readRDS("/home/labgroups/ccgg/downloads/pub/wu_2022/Immunoediting/data/all_mut_ccf_tpm.rds")
TMB<- table(all_mut_ccf$sample)
es_ccf$TMB<- TMB[es_ccf$sample]

# Get NeoAg proportions
neo_t<- table(all_mut_ccf$sample, all_mut_ccf$neo)
prop_neo<- prop.table(neo_t,1)[,"neo"]
es_ccf$prop_neo<- prop_neo[es_ccf$sample]

# Add cancer type
load("../immunoediting_2019/data/TCGA_manifest.RData")
es_ccf$cancer<- TCGA_cancer_id[substr(es_ccf$sample,1,15),1] 

# By cancer
es_cancer<- sort(tapply(es_ccf$es, es_ccf$cancer, "median", na.rm=T))
tmb_cancer<- sort(tapply(es_ccf$TMB, es_ccf$cancer, "median", na.rm=T))
neo_cancer<- sort(tapply(es_ccf$prop_neo, es_ccf$cancer, "median", na.rm=T))

cancers<- names(es_cancer)
cancer_df<- data.frame(
  cancer=cancers,
  es=es_cancer[cancers],
  tmb=tmb_cancer[cancers],
  neo=neo_cancer[cancers]
)

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

# Demonstrate that filtering mainly affects samples with low MB
#################################################################

samples<- unique(all_mut_ccf$sample)
n_samples_preFilter<- table(TCGA_cancer_id[substr(samples,1,15),1])
n_samples_postFilter<- table(es_ccf$cancer)
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
  stat_cor(method = "spearman", label.x = 200, label.y =  floor(10*min(cancer_df$pct_filter, na.rm=T))/10, hjust=1, vjust=0, size=2)

# Combine plots
p_cancer_tmb<- plot_grid(
  p_tmb_es, p_tmb_neo, p_tmb_filter,
  ncol=3,
  labels="AUTO") 

# Demonstrate that filtering low TMB cancers affects proportions of neoAg mutations
##################################################################################################

# Calculate Probablity no NeoAgs/subclonal ~ TMB
neo_t<- table(all_mut_ccf$neo)
Pneo<- prop.table(neo_t)["neo"] # 0.093

sim_t<- as.data.frame(matrix(NA, 100, 2, dimnames = list(1:100, c("TMB","neo"))))
for(i in 1:100){
  sim_t[i,"TMB"]<- i
  sim_t[i,"neo"]<- dbinom(0, size=i, prob=Pneo)
}
sim_t$col<- "black"
sim_t$col[sim_t$TMB%in%c(5,10,20,50)]<- "red"
sim_t$label<- NA
sim_t$label[sim_t$TMB%in%c(5,10,20,50)]<- c("TMB = 5", "TMB = 10", "TMB = 20", "TMB = 50")

# Plot Probablity no NeoAgs ~ TMB
sim_t<- sim_t%>%arrange(col)
p_sim_neo<- ggplot(sim_t, aes(x=TMB, y=neo, label=label)) + 
  geom_point(color=sim_t$col)+
  geom_line()+
  theme_prism()+
  ylab("Probability (number of neoantigenic mutations = 0)") +
  theme(
    # legend.position="none",
    axis.line = element_line(size=0.4),
    axis.ticks = element_line(size=0.4),
    text = element_text(size = 7)
  ) +
  ggrepel::geom_text_repel(size=2, col="red")

# Plot binomial distributions for neoantigens
for(TMB in c(5,10,20,50)){
  
  success=(0:TMB)
  binom_df<- data.frame(
    success=success,
    dbinom=dbinom(success, size=TMB, prob=Pneo)
  )
  
  p<- ggplot(binom_df,aes(x=success,xend=success,y=0,yend=dbinom,color=(success==0))) +
    geom_segment(size=1) +
    scale_color_manual(values=c("black","red")) +
    scale_x_continuous(breaks=seq(0,10,1), limits=c(0,10)) +
    # scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,ceiling(10*max(binom_df$dbinom))/10)) +
    # ggtitle(paste0("TMB = ", TMB)) +
    annotate("text", x = 10, y = ceiling(10*max(binom_df$dbinom))/10, label = paste0("TMB = ", TMB), hjust=1, vjust=1, size=2) +
    xlab("Number of neoantigenic mutations") +
    ylab("Probability") +
    theme_prism()+
    theme(
      legend.position="none",
      axis.line = element_line(size=0.4),
      axis.ticks = element_line(size=0.4),
      text = element_text(size = 7)
    )
  
  assign(paste0("p_binom_neo",TMB),p)
}

# Merge plots
p_binom_neo<- plot_grid(
  p_binom_neo5, p_binom_neo10, p_binom_neo20, p_binom_neo50,
  ncol = 2) 

p_binom_neo<- plot_grid(
  p_binom_neo, p_sim_neo,
  ncol=2,
  rel_widths = c(2,1)
  ) + theme(plot.margin = unit(c(1,0,0,1), "lines"))

p<- plot_grid(
  p_cancer_tmb, p_binom_neo,
  ncol=1,
  labels=c(NA,"D")
)

# Save
######
saveRDS(p, "results/data/p_ccf.rds")


