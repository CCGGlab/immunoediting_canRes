#############################################################################
# ES_RNA_Wu permutation 
#############################################################################

library(tidyverse)
library(parallel)
source("scripts/functions/cales.R")
source("scripts/functions/cales_t.R")

# Exclude samples with lacking expression data
mut_dt<- readRDS(file = "temp/mut_df_wu.rds")
samples<- unique(mut_dt$sample)

# Permutation
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
  saveRDS(res, paste0("temp/wu/perm_",s,".rds"))
}
mclapply(1:length(samples), "perm_es", mc.cores = 50)

# Merge permutation data in matrix 
perm_matrix<- matrix(NA,2000,length(samples))
colnames(perm_matrix)<- samples
for(f in list.files("temp/wu/", full.names = T)){
  s<- gsub(".*perm_|\\.rds","",f)
  perm_matrix[,s]<- as.numeric(readRDS(f))
}
saveRDS(perm_matrix, file = "temp/perm_matrix_wu.rds")
