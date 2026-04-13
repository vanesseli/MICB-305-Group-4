library(tidyverse)
library(phyloseq)
library(readxl)
library(vegan)
library(ggpubr)
library(microbiome)
library(indicspecies)
library(writexl)
library(ANCOMBC)
library(dplyr)
library(readr)

#loading data
metadata = read.csv('filtering/metadata_filtered.csv')
#61

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_id

#indicator species analysis
ps = readRDS("project_analyses/my_phyloseq_object.rds")|>
  tax_glom("Genus")

ps_phylum = tax_glom(ps, 'Genus')
ps_relab = transform(ps_phylum, 'compositional')
ps_filt = filter_taxa(ps_relab, function(x) mean(x) > 0.001, TRUE) #keep that is true for
otu_table = data.frame(otu_table(ps_filt))

set.seed(400)
indval = multipatt(t(otu_table), cluster = ps_filt@sam_data$antidepressant_on_off, control = how(nperm = 999))

summary(indval, indvalcomp = TRUE)
indval_table= as.data.frame(indval$sign)

#matching actual names to the codes
# extract the taxonomy table as a data frame
tax_info <- as.data.frame(tax_table(ps_filt))

# match the ISA results with the taxonomy names
indval_table$Genus <- tax_info[rownames(indval_table), "Genus"]

# move Genus to the first column so i can read
indval_table <- indval_table[, c("Genus", setdiff(names(indval_table), "Genus"))]

View(indval_table)
write_csv(indval_table, "indicator_species/indval_table.csv")



#what are the taxonomic ranks of CAG 
tax_table(ps)[grep("CAG-873", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of g__Catenibacterium
tax_table(ps)[grep("Catenibacterium", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of  g__Anaerovibrio
tax_table(ps)[grep("Anaerovibrio", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of g__Clostridia_UCG-014 
tax_table(ps)[grep("Clostridia_UCG-014", tax_table(ps)[, "Genus"]), ]

#what are the taxonomic ranks of   g__Fusobacterium
tax_table(ps)[grep("Fusobacterium", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of  g__Rikenellaceae_RC9_gut_group
tax_table(ps)[grep("Rikenellaceae_RC9_gut_group", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of  g__Muribaculaceae
tax_table(ps)[grep("Muribaculaceae", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of  g__Geobacillus
tax_table(ps)[grep("Geobacillus", tax_table(ps)[, "Genus"]), ]


#plotting
phyla_to_plot = indval_table%>%
  filter(s.off ==1 & s.on==0)%>%
  rownames()

df_of_taxa = prune_taxa(phyla_to_plot, ps_filt)%>%
  psmelt()

library(ggplot2)
plot_filter <- df_of_taxa%>%
  ggplot(aes(antidepressant_on_off, Abundance, fill=antidepressant_on_off))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=0, width=0.2)+
  facet_wrap(~Genus, ncol=4, scales = 'free')

library(scales)
indicator_species_plot <- plot_filter +
  scale_y_continuous(trans = pseudo_log_trans())+ scale_fill_manual(values = c("off" = "#4876FF","on"  = "#FF82AB" )) 

indicator_species_plot

ggsave("indicator_species/indicator_species.png", plot = indicator_species_plot, width = 10, height = 10, limitsize = FALSE)
