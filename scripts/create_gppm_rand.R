library("GenomicRanges") # version 1.28.6

# Required packages to work with this old GPPM object:
# bioconductor-genomicranges-1.28.6
# bioconductor-biostrings-2.44.2

# Load from GPPM
# These data are available at https://zenodo.org/record/2621365/files/GPPM_inclHLAAlleles.rds
GPPM <- readRDS("../immunoediting_2019/data/GPPM_inclHLAAlleles.rds")
idx_miss <- which(GPPM$variant=="nonsynonymous SNV") # Only missense
GPPM <- GPPM[idx_miss]
idx_rand <- sample(1:length(GPPM),100000) # 100000 random mutations
GPPM <- GPPM[idx_rand]

# Displaying the GenomicRanges object or converting it to a data frame using
# as.data.frame failed (even using the same version of GenomicRanges it was created with)
# Manually convert it to a dataframe:
df <-
  data.frame(
    seqnames = GPPM@pos_runs@seqnames,
    GPPM@pos_runs@ranges,
    strand = GPPM@pos_runs@strand,
    GPPM@elementMetadata
  )

# Check:
# The column "gene" was part of the metadata (elementMetadata)
# as the chromosome (seqnames), start and end position (ranges) and
# gene names seem to correspond, the converted data frame looks fine

saveRDS(df, "data/GPPM_rand_df.rds")
