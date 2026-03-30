library(tidyverse)
library(ANCOMBC)
library(Maaslin2)
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

write.csv(statistical_table, "stat_table.csv")

#39 are significant if you ignore robust
taxa_to_plot = statistical_table |>
  filter(diff_antidepressant_on_offon==T)

taxa_to_plot |>
  ggplot(aes(taxon,lfc_antidepressant_on_offon)) +
  geom_col() +
  coord_flip()








#Maaslin2 data is all insignificant after filtering for Q-value

# otu_df <- data.frame(ps_glom@otu_table)
# otu_tss <- otu_df / rowSums(otu_df)
# 
# out = Maaslin2(input_data = otu_tss,
#                input_metadata = data.frame(ps_glom@sam_data),
#                output = 'to_delete',
#                fixed_effects = c('antidepressant_on_off'),
#                reference = 'antidepressant_on_off,off',
#                normalization = 'NONE',
#                transform = 'AST',
#                min_abundance = 0.001,
#                min_prevalence = 0.1,
#                max_significance = 0.05,
#                plot_heatmap = FALSE,
#                plot_scatter = FALSE
# )
# 
# statistical_table = out$results
# 
# taxa_to_plot = statistical_table |> 
#   filter(qval < 0.05)
# 
# #readable genus names ~~~
# tax_df <- data.frame(tax_table(ps_glom)) |>
#   rownames_to_column('feature')
# 
# taxa_to_plot <- taxa_to_plot |>
#   left_join(tax_df |> select(feature, Genus), by = 'feature') |>
#   mutate(Genus = ifelse(is.na(Genus), feature, Genus))

