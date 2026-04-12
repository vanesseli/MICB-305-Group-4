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

#loading data
metadata_edited = read_excel('filtering/depression_metadata_manuscript.xlsx')
#572

metadata_edited <- as.data.frame(metadata_edited)
rownames(metadata_edited) <- metadata_edited$sample_id

metadata_clean_edited <- metadata_edited |>
  filter(str_detect(library_name, "blank", negate = TRUE))|>
  filter(!is.na(library_name))

metadata_clean_edited <- metadata_clean_edited %>%
  mutate(antidepressant_on_off = recode(as.character(antidepressant_on_off), "0" = "Off","1" = "On")) 

#control for substance use (paper did)
metadata_clean_edited <- metadata_clean_edited|>  
  filter(current_any_substance_dx == "NO")
#427 samples

metadata_clean_edited <- metadata_clean_edited|>
  filter(bdi_group != "NA")
#364

##subsetting
#people with HIV, no HCV, depressed, use antidepressants
hiv_depressed_antidepressant <- metadata_clean_edited|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == "On")
#38 people

#people with HIV, no HCV, depressed, no antidepressants
hiv_depressed_no_anti <- metadata_clean_edited|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == "Off")
#23 people in total

View(metadata_clean_edited)
#total number of people who has HIV and depression (and no HCV)
total_people <- metadata_clean_edited|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')
#61 people in total

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

View(indval_table)

#matching actual names to the codes
# extract the taxonomy table as a data frame
tax_info <- as.data.frame(tax_table(ps_filt))

# match the ISA results with the taxonomy names
indval_table$Genus <- tax_info[rownames(indval_table), "Genus"]

# move Genus to the first column so i can read
indval_table <- indval_table[, c("Genus", setdiff(names(indval_table), "Genus"))]

#what are the taxonomic ranks of CAG 
tax_table(ps)[grep("CAG-873", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of g__Catenibacterium
tax_table(ps)[grep("Catenibacterium", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of  g__Anaerovibrio
tax_table(ps)[grep("Anaerovibrio", tax_table(ps)[, "Genus"]), ]
#what are the taxonomic ranks of g__Clostridia_UCG-014 
tax_table(ps)[grep("Clostridia_UCG-014", tax_table(ps)[, "Genus"]), ]

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

ggsave("indicator_species/indicator_species.png", plot = indicator_species_plot, width = 10, height = 10, limitsize = FALSE)
