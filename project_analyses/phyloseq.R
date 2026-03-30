#phyloseq object creation

library(tidyverse)
library(phyloseq)
library(BiocManager)

#Loading in the stuff
set.seed(123)
metadata <- read.delim("depression_metadata_full.txt", row.names = 1)
metadata <- metadata |> 
  drop_na(antidepressant_on_off, antidepressant_count, 
          bdi_group, bdi_ii, hiv_status_clean, hcv) |>
  mutate(antidepressant_on_off = recode(antidepressant_on_off,
                                        `0` = "off", `1` = "on")) |>
  filter(hcv == "NO") |>
  filter(bdi_group != 'minimal') |>
  filter(hiv_status_clean == "HIV+") |>
  filter(Assay.Type == "AMPLICON")

metadata$antidepressant_on_off <- as.factor(metadata$antidepressant_on_off)

taxonomy = read.delim('taxonomy.tsv', row.names = 1)
counts = read.delim('feature-table.txt', skip=1, row.names=1)
tree = read_tree("tree.nwk")

otu_table <- otu_table(counts, taxa_are_rows = TRUE)

sample_data <- sample_data(metadata)

taxa_table <- taxonomy |>
  separate(col = Taxon, 
           into = c('Domain','Phylum','Class','Order',
                    'Family','Genus','Species'),
           sep=';', fill='right') |> 
  select(-Confidence)  |>
  as.matrix()

tax_table <- tax_table(taxa_table)

ps0 <- phyloseq(otu_table, sample_data, tax_table, tree)

saveRDS(ps0, "my_phyloseq_object.rds")
