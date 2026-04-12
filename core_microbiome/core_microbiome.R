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

##core microbiome
ps = readRDS("project_analyses/my_phyloseq_object.rds")|>
  tax_glom("Genus")
#subset phyloseq
on_anti = subset_samples(ps, total_people$antidepressant_on_off == 'On')
off_anti = subset_samples(ps, total_people$antidepressant_on_off == 'Off')

#_rare_relab_genus
#core members
core_on_anti = core_members(on_anti, detection = 0.005, prevalence = 0.3)
core_off_anti = core_members(off_anti, detection = 0.005, prevalence = 0.3)

library(ggVennDiagram)

venn <- ggVennDiagram(list(core_on_anti, core_off_anti), set_size = 5, category.names = c('On', 'Off'))
venn <- ggVennDiagram(
  list(core_on_anti, core_off_anti),
  set_size = 5,
  category.names = c('On', 'Off')
) +
  scale_fill_gradient(low = "#4876FF", high = "#FF82AB")
venn                                                                                                                   
ggsave("indicator_species/core_microbiome_plot.png", plot = venn, width = 10, height = 10, limitsize = FALSE)
