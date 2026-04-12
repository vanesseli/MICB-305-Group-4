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

#loading the data
metadata = read_excel('filtering/depression_metadata_manuscript.xlsx')
#572

metadata <- metadata |>
  filter(str_detect(library_name, "blank", negate = TRUE))|>
  filter(!is.na(library_name))

metadata <- metadata %>%
  mutate(antidepressant_on_off = recode(as.character(antidepressant_on_off), "0" = "Off","1" = "On")) 

metadata <- metadata|>
  filter(str_detect(name_sample, "blank", negate = TRUE))|>
  filter(!is.na(name_sample))

#control for substance use (paper did)
metadata <- metadata|>  
  filter(current_any_substance_dx == "NO")
#427 samples

#checking for duplicates (all of them should say false)
any(duplicated(metadata$name_sample))
any(duplicated(metadata$sample_name))
any(duplicated(metadata$BioSample))
any(duplicated(metadata$Sample_name))
any(duplicated(metadata$host_subject_id))

#if one of them return with ‘TRUE’, check which samples were duplicated with this code
#metadata_clean_edited$host_subject_id[duplicated(metadata_clean_edited$host_subject_id)]

sum(metadata$bdi_group == "NA")
#63 samples have NA for bdi group

metadata <- metadata|>
  filter(bdi_group != "NA")
#364

#want to see how many NAs there are for antidepressant column (6)
metadata <- metadata |>
  filter(antidepressant_on_off != "NA")
#359

#people with HIV, no HCV, depressed, use antidepressants
hiv_depressed_antidepressant <- metadata|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == "On")
#38 people

#people with HIV, no HCV, depressed, no antidepressants
hiv_depressed_no_anti <- metadata|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == "Off")
#23 people in total

#total number of people who has HIV and depression (and no HCV)
total_people <- metadata|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')
#61 people in total

write.csv(total_people, "filtering/metadata_filtered.csv")
