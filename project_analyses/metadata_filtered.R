metadata = read.delim('depression_metadata_full.txt', row.names = 1)
#1415 samples

metadata <- metadata|>
  filter(Assay.Type == "AMPLICON")
#1031
#amplicon because we are interested in 16s, they also did whole genome seq (WGS)

#in the prior_visit_exam_date column some are labelled ‘1’, which i believe means it is their first time coming. So instead I replace ‘1’ with NA (indicating its their first time here)
metadata_clean <- metadata |>
  dplyr::mutate(
    prior_visit_exam_date = stringr::str_replace_all(prior_visit_exam_date, "1", NA_character_))

#filtering for all the NA values (because we want to keep the first sample if they came multiple times, this is what the paper did
metadata_clean <- metadata_clean|>
  filter(is.na(prior_visit_exam_date))

#in sample name, some are ‘blank samples’ and not associated with ‘human id’
metadata_clean <- metadata_clean|>
  filter(str_detect(Sample_name, "qiita_sid", negate = TRUE))|>
  filter(!is.na(Sample_name))
#637

#in name_sample, some are ‘blank samples’ and not associated with ‘human id’
metadata_clean <- metadata_clean|>
  filter(str_detect(name_sample, "blank", negate = TRUE))|>
  filter(!is.na(name_sample))
#592

#in extractionkit_lot..exp, some are ‘blank samples’ and not associated with ‘human id’
metadata_clean <- metadata_clean|>
  filter(str_detect(extractionkit_lot..exp., "KOLStudy1.0", negate = TRUE))|>
  filter(!is.na(extractionkit_lot..exp.))
#573

#checking for duplicates in the name_sample column
any(duplicated(metadata_clean$name_sample))
metadata_clean$name_sample[duplicated(metadata_clean$name_sample)]
#6 duplicates

#removed the 6 duplicates
metadata_clean <- metadata_clean|>
  distinct(name_sample, .keep_all = TRUE)|>
  filter(!is.na(name_sample))
#567

#checking for duplicates in different columns
metadata_clean <- metadata_clean|>
  distinct(Sample.Name, .keep_all = TRUE)|>
  filter(!is.na(Sample.Name))
#no duplicates in this column, 567

metadata_clean <- metadata_clean|>
  distinct(Sample_name, .keep_all = TRUE)|>
  filter(!is.na(Sample_name))
#no duplicates in this column, 567

metadata_clean <- metadata_clean|>
  distinct(orig_name..exp., .keep_all = TRUE)|>
  filter(!is.na(orig_name..exp.))
#no duplicates in this column, 567

metadata_clean <- metadata_clean|>
  distinct(BioSample, .keep_all = TRUE)|>
  filter(!is.na(BioSample))
#no duplicates in this column, 567

metadata_clean <- metadata_clean|>
  distinct(Experiment, .keep_all = TRUE)|>
  filter(!is.na(Experiment))
#no duplicates in this column, 567

#filtering out duplicate host_subject_id
metadata_clean <- metadata_clean|>
  distinct(host_subject_id, .keep_all = TRUE)|>
  filter(!is.na(host_subject_id))
#402 samples

#there were a couple ‘samples’ that were, again, not associated with ‘humans’ so we cleaned them out
metadata_clean <- metadata_clean|>  
  filter(str_detect(host_subject_id, "Sample", negate = TRUE))|>
  filter(str_detect(host_subject_id, "T0", negate = TRUE))
#398

metadata_clean$antidepressant_on_off = as.factor(metadata_clean$antidepressant_on_off)

#remove all the NA values from relevant columns
metadata_clean <- metadata_clean%>%
  filter(!is.na(bdi_group))%>%
  filter(!is.na(antidepressant_on_off))%>%
  filter(!is.na(hiv_status_clean))%>%
  filter(!is.na(hcv))

metadata_clean <- metadata_clean %>%
  mutate(antidepressant_on_off = recode(as.character(antidepressant_on_off), "0" = "Off","1" = "On")) 

metadata_clean <- metadata_clean%>%
  filter(!is.na(antidepressant_on_off))

#people with HIV, no HCV, depressed, use antidepressants
hiv_depressed_antidepressant <- metadata_clean|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == 'On')
#34 people in total

#people with HIV, no HCV, depressed, no antidepressants
hiv_depressed_no_anti <- metadata_clean|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == 'Off')
#23 people in total

#total number of people who has HIV and depression (and no HCV)
metadata_filtered <- metadata_clean|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')
#57 people in total

write.csv(metadata_filtered, "metadata_filtered.csv")
