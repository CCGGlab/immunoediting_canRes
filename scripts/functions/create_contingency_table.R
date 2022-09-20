create_contingency_table <- function(gs_genes, category_genes, all_genes) {
  table(factor(all_genes %in% category_genes, c(FALSE, TRUE)),
        factor(all_genes %in% gs_genes, c(FALSE, TRUE)))
}