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

#These are the two identified through Maaslin2
#non-significant P-value
notable <- statistical_table |>
  filter(taxon == 'dc64f4ef4d4f99019072896947b19e4a' | taxon == '4499b0901e611126af0d04752c158976')

#this one is oribacterium (shown to have a strong lfc value)
tax_table(ps_glom)['dc64f4ef4d4f99019072896947b19e4a', ]
#this one is lachnospiraceae, a serotype of it I think
tax_table(ps_glom)['4499b0901e611126af0d04752c158976', ]

#39 are significant if you ignore robust
taxa_to_plot = statistical_table |>
  filter(diff_antidepressant_on_offon==T)

tax_df <- data.frame(tax_table(ps_glom)) |>
  rownames_to_column('taxon')

taxa_to_plot <- taxa_to_plot |>
  left_join(tax_df |> select(taxon, Genus), by = 'taxon') |>
  mutate(Genus = ifelse(is.na(Genus), taxon, Genus))

#order it from smallest to largest LFC
taxa_to_plot |>
  mutate(Genus = reorder(Genus, lfc_antidepressant_on_offon)) |>
  ggplot(aes(x = Genus, y = lfc_antidepressant_on_offon)) +
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

