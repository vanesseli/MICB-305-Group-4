# MICB-305-Group-4
General Information
This repository contains all codes required to run the analyses in the

Team members: Vanesse Li, Esther Xu, Ray Lou, Sarah Lai, Katelyn Wang

##Research Question: Does antidepressant use lead to additional changes in gut microbiome composition, diversity, and predicted functional pathways in individuals with HIV and major depressive disorder (HIV+MDD+)?

Methods:
- alpha diversity analysis
- beta diversity analysis
- indicator species analysis
- differential abundance analysis
- PICRUSt2 predictive functional analysis
- random forest

Details on the dataset:
- "depression_metadata_with_rationale" includes the full dataset & rationale on why specific samples were not included (ex: blanks, patient duplicates, errors, etc)
- "depression_metadata_manuscript" excludes the duplicates and blanks. patient without information on HIV status, HCV status, antidepressant status, and depression status are also not included in here. To get to 'metadata_filtered' (our cohort of interest), please run the 'filtering.R' Rscript
- "metadata_filtered" contains the actual 61 patient data we are working with
  
