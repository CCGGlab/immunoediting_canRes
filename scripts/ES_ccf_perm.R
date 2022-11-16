#############################################################################
# ES_RNA permutation 
#############################################################################


library(tidyverse)
library(parallel)
source("scripts/functions/cal_es_new_test.R")

# Load
mut_dt<- readRDS(file = "temp/mut_df_ccf.rds")
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
    # cat(j, " ")
    mut_dt_perm$neo<- sample(mut_dt_perm$neo)
    # Calculate ES
    ES_ccf_perm<- unlist(cal_es_new_test( mut_dt_perm))
    res<- append(res, ES_ccf_perm)
  }
  # Save 
  saveRDS(res, paste0("temp/ccf/perm_ccf",s,".rds"))
}
mclapply(1:length(samples), "perm_es", mc.cores = 50)

# Merge permutation data in matrix 
perm_matrix<- matrix(NA,2000,length(samples))
colnames(perm_matrix)<- samples
for(f in list.files("temp/ccf/", full.names = T)){
  s<- gsub(".*ccf|\\.rds","",f)
  perm_matrix[,s]<- as.numeric(readRDS(f))
}
saveRDS(perm_matrix, file = "temp/perm_matrix_ccf.rds")
