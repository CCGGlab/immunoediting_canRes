library(tidyverse)

# ==== Load auxiliary functions ====

source("scripts/functions/read_go.R")
source("scripts/functions/read_gtex.R")

# ==== Read data ====

rand <- readRDS("data/GPPM_rand_df.rds")

# GTex expression data
expression <- read_gtex()

# Gene Ontology data
go <- deframe(read_go())
go_molecular_functions <-
  list(GCR = "GOMF_G_PROTEIN_COUPLED_RECEPTOR_ACTIVITY",
       OLFR = "GOMF_OLFACTORY_RECEPTOR_ACTIVITY") %>%
  # Get list of genes for these gene set names
  map( ~ go[[.]])

# Define hydrophobic amino acids
hydrophobic <- c("A", "F", "G", "I", "L", "M", "P", "V", "W")

# ==== Analysis ====

mutations <-
  transmute(
    rand,
    gene,
    # Kd < 500: binding
    mhc1_bind = HLA_aff_mean < 500,
    # Calculate proportion of hydrophobic amino acids
    hydrophobic_proportion = map_dbl(nonaPep_wt, function(nonapeps) {
      # Calculate which proportion of amino acids
      # in all possible 9-mers is hydrophobic
      # Convert nonapeptide to list of amino acids
      amino_acids <-
        nonapeps %>%
        paste0(collapse = "") %>%
        str_split("") %>%
        unlist
      num_hp <- sum(amino_acids %in% hydrophobic)
      num_all <- length(amino_acids)
      num_hp / num_all
    })
  ) %>%
  # Gene occurs in molecular function
  # -> Boolean column per molecular function
  mutate(map_dfc(go_molecular_functions, ~ gene %in% .))

# Combine mutation and expression value
mut_expr <- left_join(mutations, expression, by="gene")

# Threshold expression values into not expressed and expressed
mut_expr_thesholded <- mut_expr %>%
  # For 1309 mutations no expression value was found
  filter(!is.na(expression)) %>%
  mutate(is_expressed = expression > 0)

# Calculate p values
fisher_bind_expr <- fisher.test(table(mut_expr_thesholded$mhc1_bind, mut_expr_thesholded$is_expressed))
fisher_bind_olfr <- fisher.test(table(mut_expr_thesholded$mhc1_bind, mut_expr_thesholded$OLFR))
fisher_bind_gcr <- fisher.test(table(mut_expr_thesholded$mhc1_bind, mut_expr_thesholded$GCR))

# Save to RDS
saveRDS(mut_expr_thesholded, "data/gppm_mut_expr_thresholded.rds")
saveRDS(list(expr=fisher_bind_expr, olfr=fisher_bind_olfr, gcr=fisher_bind_gcr), "data/gppm_nonexp_pvalues.rds")

# ==== GSEA barplots & logistic regression==== 
######################################################

# Load RDS
mut_expr_thesholded <- readRDS("data/gppm_mut_expr_thresholded.rds")
fisher_test_pvalues <- readRDS("data/gppm_nonexp_pvalues.rds")

# Make plots and store to list
plts <- list(
  expression = mut_expr_thesholded %>%
    # Calculate proportion of predicted neoantigens for the groups
    # gene is expressed: yes / no"
    group_by(is_expressed) %>%
    summarise(predicted_neoantigens = mean(mhc1_bind)) %>%
    ggplot(aes(x = is_expressed, y = predicted_neoantigens)) +
    geom_col() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.05), limits = c(0,0.35)) +
    # Add p-values (based on Fisher test) to the plot
    annotate(
      "text",
      x = 1.5,
      y = 0.31,
      hjust = 0.5,
      vjust=1,
      size=2,
      label = str_glue("p = {format(fisher_test_pvalues$expr$p.value, digits=3)}")
    ) +
    labs(x = "expressed"),
  
  olfr = mut_expr_thesholded %>%
    # Calculate proportion of predicted neoantigens for the groups
    # "gene is a olfactory receptor: yes / no"
    group_by(OLFR) %>%
    summarise(predicted_neoantigens = mean(mhc1_bind)) %>%
    ggplot(aes(x = OLFR, y = predicted_neoantigens)) +
    geom_col() +
    scale_y_continuous(limits = c(0,0.55)) +
    # Add p-values (based on Fisher test) to the plot
    annotate(
      # geom_text(position = position_nudge(x = -1)),
      "text",
      x = 1.5,
      y = 0.55,
      hjust = 0.5,
      vjust=1,
      size=2,
      label = str_glue("p = {format(fisher_test_pvalues$olfr$p.value, digits=3)}")
    ),
  
  gcr = mut_expr_thesholded %>%
    # Calculate proportion of predicted neoantigens for the groups
    # "gene is a G coupled receptor: yes / no"
    group_by(GCR) %>%
    summarise(predicted_neoantigens = mean(mhc1_bind)) %>%
    ggplot(aes(x = GCR, y = predicted_neoantigens)) +
    geom_col() +
    scale_y_continuous(limits = c(0,0.45)) +
    # Add p-values (based on Fisher test) to the plot
    annotate(
      # geom_text(position = position_nudge(x = -1)),
      "text",
      x = 1.5,
      y = 0.48,
      hjust = 0.5,
      vjust=1,
      size=2,
      label = str_glue("p = {format(fisher_test_pvalues$gcr$p.value, digits=3)}")
    ),
  
  hydrophobic = ggplot(mut_expr_thesholded, aes(x = hydrophobic_proportion, y = as.numeric(mhc1_bind))) +
    # Plot logistic regression line (glm model with logit linker function)
    stat_smooth(
      method = "glm",
      se = FALSE,
      method.args = list(family = binomial)
    ) +
    labs(x = "hydrophobic proportion")
) %>%
  # Apply theme and axis labels that are common for all plots
  map(function(x) {
    x + ggpubr::theme_pubclean() +
      theme(axis.title = element_text(size = 10)) +
      labs(y = "predicted neoantigens")
  })

# Save resulting plots to RDS
saveRDS(plts, "results/data/GSEA_barplots_ggplot.rds")
