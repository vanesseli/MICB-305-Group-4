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
metadata_edited = read_excel('filtering/depression_metadata_manuscript.xlsx')
#572

metadata_clean_edited <- metadata_edited |>
  filter(str_detect(library_name, "blank", negate = TRUE))|>
  filter(!is.na(library_name))

metadata_clean_edited <- metadata_clean_edited %>%
  mutate(antidepressant_on_off = recode(as.character(antidepressant_on_off), "0" = "Off","1" = "On")) 

metadata_clean_edited <- metadata_clean_edited|>
  filter(str_detect(name_sample, "blank", negate = TRUE))|>
  filter(!is.na(name_sample))

#control for substance use (paper did)
metadata_clean_edited <- metadata_clean_edited|>  
  filter(current_any_substance_dx == "NO")

#427 samples

#checking for duplicates (all of them should say false)
any(duplicated(metadata_clean_edited$name_sample))
any(duplicated(metadata_clean_edited$sample_name))
any(duplicated(metadata_clean_edited$BioSample))
any(duplicated(metadata_clean_edited$Sample_name))
any(duplicated(metadata_clean_edited$host_subject_id))

#if one of them return with ‘TRUE’, check which samples were duplicated with this code
#metadata_clean_edited$host_subject_id[duplicated(metadata_clean_edited$host_subject_id)]

#only select the columns I care about since the dataset is so wide
metadata_clean_small <- metadata_clean_edited |>
  select(hiv_status_clean, antidepressant_on_off, hcv, bdi_group, host_subject_id, current_any_substance_dx)
#427

sum(metadata_clean_small$bdi_group == "NA")
#63 samples have NA for bdi group

metadata_clean_small <- metadata_clean_small|>
  filter(bdi_group != "NA")
#364

#want to see how many NAs there are for antidepressant column (6)
metadata_clean_small <- metadata_clean_small|>
  filter(antidepressant_on_off != "NA")
#359

metadata_clean_small_1 <- metadata_clean_small|>
  filter(hcv == "NO")|>
  filter(hiv_status_clean=="HIV+")|>
  filter(bdi_group!="minimal")|>
  filter(antidepressant_on_off == "On")

#people with HIV, no HCV, depressed, use antidepressants
hiv_depressed_antidepressant <- metadata_clean_small|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == "On")
#38 people

#people with HIV, no HCV, depressed, no antidepressants
hiv_depressed_no_anti <- metadata_clean_small|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == "Off")
#23 people in total

#total number of people who has HIV and depression (and no HCV)
total_people <- metadata_clean_small|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')
#61 people in total
