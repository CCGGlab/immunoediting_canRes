read_gtex <- function() {
  read_tsv(
    "downloads/GTex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
    skip = 2
  ) %>%
    dplyr::rename(ensembl = Name, gene = Description) %>%
    pivot_longer(-c(ensembl, gene),
                 names_to = "tissue",
                 values_to = "expression") %>%
    # Median expression per gene (over all tissue types)
    group_by(gene) %>%
    summarise(expression = median(expression, na.rm = T))
}
