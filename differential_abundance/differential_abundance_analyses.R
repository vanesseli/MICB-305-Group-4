library(tidyverse)
library(ANCOMBC)
library(Maaslin2)
library(phyloseq)
library(ggplot2)
library(readxl)

set.seed(421)

ps_glom = readRDS("my_phyloseq_object.rds") |>
  tax_glom("Genus")

#Retain samples with reads more than 1000
ps_filt <- prune_samples(sample_sums(ps_glom) > 1000, ps_glom)









#ANCOMBC2 analysis
out_ancom = ancombc2(data = ps_filt,
                     fix_formula = 'antidepressant_on_off',
                     p_adj_method = 'BH',
                     prv_cut = 0.1)

ancom_statistical_table = out_ancom$res

#41 (37) are significant if you ignore robust
taxa_to_plot = ancom_statistical_table |>
  filter(diff_antidepressant_on_offon==T)
count(taxa_to_plot)

#add a column for genus names
tax_df <- data.frame(tax_table(ps_filt)) |>
  rownames_to_column('taxon')
taxa_to_plot <- taxa_to_plot |>
  left_join(tax_df |> select(taxon, Genus), by = 'taxon') |>
  mutate(Genus = ifelse(is.na(Genus), taxon, Genus))

#order it from smallest to largest LFC
lfc_ancom <-
  mutate(taxa_to_plot, Genus = reorder(Genus, lfc_antidepressant_on_offon)) |>
  ggplot(aes(x = Genus, y = lfc_antidepressant_on_offon, fill = "#4876FF")) +
  geom_col() +
  coord_flip()
lfc_ancom
ggsave("lfc_ancom.png", width = 5, height = 5, dpi = 300)

#Volcano updated colours
library(ggrepel)
volcano_df <- ancom_statistical_table |>
  mutate(
    log_q = -log10(q_antidepressant_on_offon),
    significant = diff_antidepressant_on_offon,
    direction = case_when(
      diff_antidepressant_on_offon & lfc_antidepressant_on_offon > 0 ~ "Up (On)",
      diff_antidepressant_on_offon & lfc_antidepressant_on_offon < 0 ~ "Down (On)",
      TRUE ~ "NS"
    )
  ) |>
  left_join(tax_df |> select(taxon, Genus), by = "taxon") |>
  mutate(
    Genus = ifelse(is.na(Genus), taxon, Genus),
    label = NA
  )
updated_colour <- ggplot(volcano_df, aes(x = lfc_antidepressant_on_offon, y = log_q)) +
  geom_point(
    data = volcano_df |> filter(direction == "NS"),
    color = "grey70", alpha = 0.3, size = 1.2
  ) +
  geom_point(
    data = volcano_df |> filter(direction != "NS"),
    aes(color = direction),
    size = 2.5, alpha = 0.85
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  scale_color_manual(
    values = c("Up (On)" = "#FF82AB", "Down (On)" = "#4876FF"),
    name = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Log Fold Change (On / Off)",
    y = expression(-log[10](q-value)),
    title = "ANCOM-BC2 Volcano Plot",
    subtitle = "Dashed line: q = 0.05 | Pink/Blue = +/- LFC"
  )
updated_colour

ggsave("updated_colour_data.png", width = 5, height = 5, dpi = 300)










#Maaslin2 Analysis
out_maaslin = Maaslin2(input_data = data.frame(ps_filt@otu_table),
                       input_metadata = data.frame(ps_filt@sam_data),
                       output = 'to_delete',
                       fixed_effects = c('antidepressant_on_off'),
                       reference = 'antidepressant_on_off,off',
                       normalization = 'TSS',
                       transform = 'AST',
                       min_abundance = 0.001,
                       min_prevalence = 0.1,
                       max_significance = 0.05,
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE
)

maaslin_statistical_table = out_maaslin$results

#no significant taxa
maaslin_taxa_to_plot = maaslin_statistical_table |>
  filter(qval < 0.05)

#renaming column name feature to taxon for the next step of adding genus names
maaslin_statistical_table <- maaslin_statistical_table |>
  rename(taxon = feature)

#some features have an X in the beginning, which is making the left join not work for some features
maaslin_statistical_table <- maaslin_statistical_table |>
  mutate(taxon = sub("^X", "", taxon))

#joining the names to the table
maaslin_taxa_to_plot <- maaslin_statistical_table |>
  left_join(tax_df|> select(taxon, Genus), by = "taxon") |>
  mutate(Genus = ifelse(is.na(Genus), taxon, Genus))


#Top 2 taxa with notable p-values noted here
#this is genus Holdemania, 0.2309500, but very little distribution (22/38)
tax_table(ps_glom)['d348b3e42b9cd325d1f0088462aedcf3', ]
#this one is lachnospiraceae, a serotype of it I think, 0.2309500, (9/38)
tax_table(ps_glom)['4499b0901e611126af0d04752c158976', ]










#Analysis

#Maaslin2 ANCOMBC overlap, only one overlap
#non-significant P-value
notable_ancom <- taxa_to_plot |>
  filter(taxon == 'd348b3e42b9cd325d1f0088462aedcf3' | taxon == '4499b0901e611126af0d04752c158976')
#g__Lachnospiraceae_NK3A20_group is now overlapping

#1 Indicator_species ANCOMBC overlap
indicator_overlap_ANCOMBC <- taxa_to_plot |>
  filter(Genus == " g__CAG-873" | Genus == " g__Catenibacterium"
         | Genus == " g__Anaerovibrio" | Genus == " g__Clostridia_UCG-014"
         | Genus == " g__Fusobacterium")
#g__CAG-873 aligns

#5 non-significant Indicator_species Maaslin2 overlap
indicator_overlap_Maaslin2 <- maaslin_taxa_to_plot |>
  filter(Genus == " g__CAG-873" | Genus == " g__Catenibacterium"
         | Genus == " g__Anaerovibrio" | Genus == " g__Clostridia_UCG-014"
         | Genus == " g__Fusobacterium")
#The taxa here that overlaps is Fusobacterium, it has a q-value of over 0.5