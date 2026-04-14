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
metadata = read.csv('filtering/metadata_filtered.csv')
#61

metadata<- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_id

metadata <- metadata%>%
  mutate(antidepressant_on_off = recode(as.character(antidepressant_on_off), "0" = "Off","1" = "On")) 

##core microbiome
ps = readRDS("project_analyses/my_phyloseq_object.rds")|>
  tax_glom("Genus")
#subset phyloseq
on_anti = subset_samples(ps, metadata$antidepressant_on_off == 'On')
off_anti = subset_samples(ps, metadata$antidepressant_on_off == 'Off')

#_rare_relab_genus
#core members
core_on_anti = core_members(on_anti, detection = 0.005, prevalence = 0.5)
core_off_anti = core_members(off_anti, detection = 0.005, prevalence = 0.5)

library(ggVennDiagram)

venn <- ggVennDiagram(list(core_on_anti, core_off_anti), set_size = 5, category.names = c('On', 'Off'))
venn <- ggVennDiagram(
  list(core_on_anti, core_off_anti),
  set_size = 5,
  category.names = c('On', 'Off')
) +
  scale_fill_gradient(low = "#4876FF", high = "#FF82AB")
venn                                                                                                                   
ggsave("indicator_species/core_microbiome_plot_0.5.png", plot = venn, width = 10, height = 10, limitsize = FALSE)
