Analysis to demonstrate that the ES RNA analysis reported in the manuscript of [Wu et al., 2022](https://aacrjournals.org/cancerres/article/82/12/2226/699353/Quantification-of-Neoantigen-Mediated) is unrelated to immunoediting.

# Simulated mutation data

Data derived from https://zenodo.org/record/2621365/files/TCGA_maf_sim.rds.
For details, see [Van den Eynden et al., 2019](https://www.nature.com/articles/s41588-019-0532-6) 

# ES RNA method

We used the *cales_t.R* method as available on https://github.com/wt12318/NeoEnrichment/

# Analysis 

```{r}
source("scripts/ES_rna_analysis")
```

results/figs/fig1.pdf

This randomly simulated mutation analysis shows that ES RNA is highly significant in both the complete (pan-cancer) dataset and for several individual cancer types. These results are strikingly similar to what Wu et al. report in fig. 2C and provides clear evidence that this result is unrelated to immunoediting (no immune selection signals are expected in randomly simulated mutation data)


