#############################################################################
# ES_RNA permutation 
#############################################################################

# 1) Load previously published simulated somatic mutation data
##############################################################

# These data are available at https://zenodo.org/record/2621365/files/TCGA_maf_sim.rds
TCGA_maf_sim<- readRDS("../immunoediting_2019/data/TCGA_maf_sim.rds")

# Only missense mutations
TCGA_maf_sim<- TCGA_maf_sim[TCGA_maf_sim$Variant_Classification=="nonsynonymous SNV",]

# 2) Calculate ES
###################

source("immunoediting_canRes/scripts/functions/cales.R")
source("immunoediting_canRes/scripts/functions/cales_t.R")

# Create mutation datatable as specified by readme in https://github.com/wt12318/NeoEnrichment/
mut_dt<- data.frame(sample=TCGA_maf_sim$Tumor_Sample_Barcode, neo="not_neo", exp=TCGA_maf_sim$mRNA)
mut_dt$neo[TCGA_maf_sim$mut_HLA_mean_aff<500]<- "neo" # Criterium to define neoantigen in Van den Eynden, 2019

# Exclude samples with lacking expression data
mut_dt<- mut_dt[!is.na(mut_dt$exp),]
samples<- unique(mut_dt$sample)

# Calculate ES RNA
library(tidyverse)
library(parallel)

# mut_dt_perm<- mut_dt
# perm_es<- function(i){
#   # Permute
#   mut_dt_perm$neo<- sample(mut_dt_perm$neo)
#   # Calculate
#   ES_rna_perm<- unlist(lapply(samples, function(x) cales_t(data = mut_dt_perm, barcode = x, calp = F, cal_type = "exp", type = "II", sample_counts = 1000)))
#   # Save 
#   saveRDS(ES_rna_perm, paste0("immunoediting_canRes/temp/perm",i,".rds"))
# }
# mclapply(1:2000, "perm_es", mc.cores = 50)

# Within each sample???
perm_es<- function(i){
  # Get sample
  s<- samples[i]
  # Get mut data for sample
  mut_dt_perm<- mut_dt[mut_dt$sample==s,]
  # Permute
  res<- NULL
  for(j in 1:2000){
    cat(j, " ")
    mut_dt_perm$neo<- sample(mut_dt_perm$neo)
    # Calculate ES
    ES_rna_perm<- unlist(cales_t(data = mut_dt_perm, barcode = s, calp = F, cal_type = "exp", type = "II", sample_counts = 1000))
    res<- append(res, ES_rna_perm)
  }
  # Save 
  saveRDS(res, paste0("immunoediting_canRes/temp/perm_",s,".rds"))
}
mclapply(1:length(samples), "perm_es", mc.cores = 50)

# Merge permumtation data in matrix 
perm_matrix<- matrix(NA,2000,length(samples))
colnames(perm_matrix)<- samples
for(f in list.files("immunoediting_canRes/temp/", full.names = T)){
  s<- gsub(".*perm_|\\.rds","",f)
  perm_matrix[,s]<- as.numeric(readRDS(f))
}
saveRDS(perm_matrix, file = "immunoediting_canRes/data/perm_matrix.rds")
