# Figure 1
###########

library(cowplot)
library(ggplot2)

# Wu 
p_wu<- readRDS("results/data/p_wu.rds")

title_wu <- ggdraw() + 
  draw_label(
    "TCGA somatic mutation data as reported in Wu et al., 2022",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5,
    size = 10
  ) 

p_wu<- plot_grid(title_wu,p_wu,ncol = 1,rel_heights = c(0.1, 1))

# Sim
p_sim<- readRDS("results/data/p_sim.rds")

title_sim <- ggdraw() + 
  draw_label(
    "Random somatic mutation data as reported in Van den Eynden et al., 2019",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5,
    size = 10
  ) 

p_sim<- plot_grid(title_sim,p_sim,ncol = 1,rel_heights = c(0.1, 1))

# GTex
p_gt<- readRDS("results/data/p_GTex.rds")

title_gt <- ggdraw() + 
  draw_label(
    "Random somatic mutations with expression data from GTEx",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5,
    size = 10
  ) 

p_gt<- plot_grid(title_gt,p_gt,ncol = 1,rel_heights = c(0.1, 1))

p<- plot_grid(
  p_wu, p_sim, p_gt,
  ncol = 1,
  labels = "AUTO")

# Footer
footer<- ggdraw() + 
  draw_label(
    "Supplementary Figure 1",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p, footer,
  ncol = 1,
  rel_heights = c(20,1)
)


ggsave("results/figs/fig1.pdf", p, width = 178, height = 265, units = "mm")
ggsave("results/figs/fig1.png", p, width = 178, height = 265, units = "mm", bg = "white")

# Figure 2
###########

p_ls<- readRDS("results/data/GSEA_barplots_ggplot.rds")

p_exp<- p_ls$expression +
  scale_y_continuous(name = "Predicted Neoantigens", labels = scales::percent) +
  scale_x_discrete(name = "Gene Expressed?", labels = c("No", "Yes")) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.margin = unit(c(1,0,0,1), "lines")
  )
  
p_olfr<- p_ls$olfr +
  scale_y_continuous(name = "Predicted Neoantigens", labels = scales::percent) +
  scale_x_discrete(name = "Olfactory Receptor?", labels = c("No", "Yes"))   +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.margin = unit(c(1,0,0,1), "lines")
  )

p_gpcr<- p_ls$gcr +
  scale_y_continuous(name = "Predicted Neoantigens", labels = scales::percent) +
  scale_x_discrete(name = "G-Protein Coupled Receptor?", labels = c("No", "Yes"))  +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.margin = unit(c(1,0,0,1), "lines")
  )

p_hp<- p_ls$hydrophobic +
  scale_y_continuous(name = "Predicted Neoantigens", labels = scales::percent) +
  scale_x_continuous(name = "Hydrophobic amino acids", labels = scales::percent)   +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.margin = unit(c(1,0,0,1), "lines")
  )

p_GSEA<- readRDS("results/data/GSEA_top5_ggplot.rds")
p_GSEA<- p_GSEA$simulated_tcga

p_aux1<- plot_grid(
  p_exp, p_gpcr, p_olfr,
  ncol = 1,
  labels = c("A","C","D"))

p_aux2<- plot_grid(
  p_GSEA+ggtitle(""), p_hp,
  ncol = 1,
  labels = c("B","E") )

p<- plot_grid(
  p_aux1, p_aux2,
  ncol = 2,
  rel_widths = c(2.2,6)
)

# Footer
footer<- ggdraw() + 
  draw_label(
    "Supplementary Figure 2",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p, NA, footer,
  ncol = 1,
  rel_heights = c(10,10,1)
)

ggsave("results/figs/fig2.pdf", p, width = 178, height = 265, units = "mm")
ggsave("results/figs/fig2.png", p, width = 178, height = 265, units = "mm")

# Figure 3
##########

p1<- readRDS("results/data/p_ccf.rds")
p2<- readRDS("results/data/p_shuffle_ccf.rds")

p<- plot_grid(
  p1, p2,
  ncol = 1
)

footer<- ggdraw() + 
  draw_label(
    "Supplementary Figure 3",
    fontface = 'bold',
    y = 0,
    hjust = 0.5,
    vjust = -1,
    size = 10
  ) 

p<- plot_grid(
  p, footer,
  ncol = 1,
  rel_heights = c(20,1)
)

ggsave("results/figs/fig3.pdf", p, width = 178, height = 265, units = "mm")
ggsave("results/figs/fig3.png", p, width = 178, height = 265, units = "mm", bg = "white")
