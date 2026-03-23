library(tidyverse)
library(ANCOMBC)
library(phyloseq)
library(ggplot2)

set.seed(421)

ps_glom = readRDS("my_phyloseq_object.rds") |>
  tax_glom("Genus")

out = ancombc2(data = ps_glom, 
               fix_formula = 'antidepressant_on_off', 
               p_adj_method = 'BH',
               prv_cut = 0.1)

statistical_table = out$res

#nothing was significant, this plot is empty
taxa_to_plot = statistical_table |>
  filter(diff_robust_antidepressant_on_offon==T)

#thought I saw some data with these parameters as true but a false robust difference, guess not though
taxa_to_plot_insignif = statistical_table |> 
  filter(diff_antidepressant_on_offon==T & passed_ss_antidepressant_on_offon==T)

statistical_table |>
  ggplot(aes(taxon,lfc_antidepressant_on_offon)) +
  geom_col() +
  coord_flip()
