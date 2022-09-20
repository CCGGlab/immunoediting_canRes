read_go <- function() {
  read_lines("downloads/mSigDB/v751/c5.go.v7.5.1.symbols.gmt") %>%
  str_split('\t') %>%
  as_tibble_col("genes") %>%
  # Extract gene set name
  mutate(gene_set = map_chr(genes, 1)) %>%
  # Remove gene set name and URL from gene list
  mutate(genes = map(genes, tail, -2)) %>%
  # Reorder columns
  select(gene_set, genes)
}
