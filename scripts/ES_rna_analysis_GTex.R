# ES_rna_analysis_GTex

library(reshape2)
library(tidyverse)
library(ggprism)

# 1) Take 100,000 random missense mutations 
###########################################

# See scripts/create_gppm_rand.R
GPPM_rand<- readRDS("data/GPPM_rand_df.rds")

# Create mut_df
mut_df<- data.frame(
  gene=GPPM_rand$gene,
  isHA=GPPM_rand$HLA_aff_mean<500
)

# 2) Add expression data from GTex
##################################

# Median expression GTex
GTex<- read.table("downloads/GTex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",sep = "\t", skip=2, header = T, check.names = F)

# merge gene expression information with mutation data
GTex$gene<- GTex$Description
GTex$Name<- NULL
GTex$Description<- NULL
mut_GTex<- merge(mut_df, GTex, by = "gene")

# Create dataframe
GTex_df<- melt(mut_GTex)

# 3) Calculate ES
###################

source("scripts/functions/cales.R")
source("scripts/functions/cales_t.R")

# Create mutation datatable as specified by readme in https://github.com/wt12318/NeoEnrichment/
mut_dt<- GTex_df
mut_dt$gene<- NULL 
mut_dt$sample<- mut_dt$variable
mut_dt$variable<- NULL 
mut_dt$neo<- "not_neo"
mut_dt$neo[mut_dt$isHA]<- "neo"
mut_dt$isHA<- NULL 
mut_dt$exp<- mut_dt$value
mut_dt$value<- NULL 

# Calculate ES RNA
samples<- unique(mut_dt$sample)
ES_rna<- lapply(samples, function(x) cales_t(data = mut_dt, barcode = x, calp = F, cal_type = "exp", type = "II", sample_counts = 1000))
names(ES_rna)<- samples

# Plot
######

# ES barplot per tissue
es_df<- data.frame(es=as.numeric(ES_rna), tissue=names(ES_rna))
es_df$tissue<- factor(es_df$tissue, es_df$tissue[levels=order(es_df$es)])
p<- ggplot(es_df, aes(tissue, es)) +  
  geom_bar(stat = "identity") +
  ylab("ES RNA") +
  xlab("") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=6),
    axis.title.y = element_text(size = 7)
  ) +
  geom_hline(yintercept=median(es_df$es),
           color = "black", linetype=2, size=0.3) +
  # geom_text(x=length(ES_rna),y=round(median(as.numeric(ES_rna), na.rm=T),digits = 3), label=paste0("Median = ",round(median(as.numeric(ES_rna), na.rm=T),digits = 3)),hjust=1,vjust=0,size=3)
  annotate(geom = "text", x=length(ES_rna),y=round(median(as.numeric(ES_rna), na.rm=T),digits = 3), label=paste0("Median = ",round(median(as.numeric(ES_rna), na.rm=T),digits = 3)),hjust=1,vjust=0,size=3)

# "Pan"?
summary(as.numeric(ES_rna))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.04877 -0.04567 -0.04282 -0.04286 -0.04074 -0.03350 
# 

# save plot
saveRDS(p, file = "results/data/p_GTex.rds")
ggsave("results/figs/fig1c.pdf", p, width = 178, height = 265/3, units = "mm")

# Save results
##############
save.image(file = "results/data/results_GTex.RData")

