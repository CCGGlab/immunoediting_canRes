library(tidyverse)

# ==== Load auxiliary functions ====

source("scripts/functions/create_contingency_table.R")
source("scripts/functions/read_go.R")
source("scripts/functions/read_gtex.R")

# ==== Read data ====

# Read and parse Gene Ontology gene sets
go <- read_go()

# Median expression GTex
gtex <- read_gtex()

# Read 3 different mutation datasets
mutations_data <- list(
  # As published by Wu, 2022
  wu = readRDS("downloads/pub/wu_2022/Immunoediting/data/all_mut_tpm_not_filter.rds"),
  # As published by Van den Eynden, 2019
  simulated_tcga = readRDS("../immunoediting_2019/data/TCGA_maf_sim.rds"),
  # Newly generated with expression data from GTeX
  gppm_rand = readRDS("data/GPPM_rand_df.rds")
)

# ==== Construct gene lists for enrichment analysis ====

# List all genes that are included in Wu's mutation database
all_genes <- list(
  wu = function(df) {
    unique(df$gene)
  },
  simulated_tcga = function(df) {
    unique(df$Hugo_Symbol)
  },
  gppm_rand = function(df) {
    unique(df$gene)
  }
)

# Genes are considered to be non-expressed if TPM <= 1 in the median patient
nonexpr_genes <- list(
  wu = function(df) {
    df %>%
      # Ensure we have a single expression value for a gene per patient
      # (This should already be the case)
      group_by(gene, sample) %>%
      summarise(tpm_exp = median(tpm_exp, na.rm = T)) %>%
      # Calculate the median expression value of the gene over the patients
      group_by(gene) %>%
      summarise(tpm_exp = median(tpm_exp, na.rm = T)) %>%
      # TPM <= 1 is considered to be not expressed
      filter(tpm_exp <= 1) %>%
      pull(gene) %>%
      unique
  },
  simulated_tcga = function(df) {
    df %>%
      # Ensure we have a single expression value for a gene per patient
      group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      summarise(mRNA = median(mRNA, na.rm = T)) %>%
      # Calculate the median expression value of the gene over the patients
      group_by(Hugo_Symbol) %>%
      summarise(mRNA = median(mRNA, na.rm = T)) %>%
      # Not expressed when mRNA == 0
      filter(mRNA == 0) %>%
      pull(Hugo_Symbol) %>%
      unique
  },
  gppm_rand = function(df) {
    # No need to summarise, this is already the median expression value
    left_join(df, gtex, by="gene") %>%
      filter(expression == 0) %>%
      pull(gene) %>%
      unique
  }
)

# ==== Perform gene set enrichment analysis ====
perform_enrichment <- function(data, nm) {
  all_genes_lst <- all_genes[[nm]](data)
  ne_genes_lst <- nonexpr_genes[[nm]](data)
  
  go %>%
    # Construct a contingency tables per gene sets in GO
    # Count for all genes under consideration in Wu study:
    # columns: gene is in GO gene set
    # rows: gene is non-expressed
    mutate(ct = map(
      genes,
      create_contingency_table,
      ne_genes_lst,
      all_genes_lst
    )) %>%
    # Perform fisher exact test and return results as columns (broom)
    mutate(map_dfr(ct, compose(broom::tidy, fisher.test))) %>%
    # Remove contingency tables from the final object
    select(-ct)
}

gsea <- mutations_data %>%
  imap(perform_enrichment)

# ==== Save results of the analysis to RDS ====
saveRDS(gsea, "results/data/GSEA_results.rds")

# ==== table GSEA top 5 ==== 
###############################

# ==== Load auxiliary functions ====

source("scripts/functions/adapt_tbl_cols.R")

# ==== Read and format results ====

# Read results of Fisher-test based enrichment analysis
gsea <- readRDS("data/GSEA_results.rds")

# Filter results that are included in the table (top 5 lowest p-value)
formatted_df <- gsea %>%
  map(function(x) {
    # Sort by p value (lowest first)
    arrange(x, p.value) %>%
      # Only display top 5 of results
      head(5) %>%
      # Give the table sensible column names
      transmute(`gene set` = gene_set,
                `p-value` = format(p.value, digits = 3, scientific=T))
  })


# ==== Plot table ====

res <- withr::with_package("ggpubr", {
  # Setup default formatting of the table
  ttheme <- ttheme(tbody.style = tbody_style(hjust = 0, x = 0.01, size=8))
  
  formatted_df %>% map(# Create the table with default formatting
    ~ ggtexttable(., theme = ttheme) %>%
      # Center align the p-value
      table_apply_fun(
        column = 3,
        fun = function(grob) {
          # Properties
          # Place the center of the text (hjust = 0.5)
          # at the center of the box (x = 0.5 npc)
          props <- list(x = unit(0.5, 'npc'), hjust = 0.5)
          # Apply properties
          reduce2(names(props), props, assign_in, .init = grob)
        }
      ))
})

# ==== Save table as RDS ====

saveRDS(res, "results/data/GSEA_top5_ggplot.rds")

