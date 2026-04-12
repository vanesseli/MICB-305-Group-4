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

metadata_edited = read_excel('depression_metadata_manuscript.xlsx')
#1415

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
  filter(str_detect(current_any_substance_dx, "YES", negate = TRUE))
#532 samples

#only select the columns I care about since the dataset is so wide
metadata_clean_small <- metadata_clean_edited |>
  select(hiv_status_clean, antidepressant_on_off, hcv, bdi_group, host_subject_id)
#532

#want to see how many NAs there are for BDI group
sum(metadata_clean_small$bdi_group == "NA")
#151 samples have NA for bdi group

metadata_clean_small <- metadata_clean_small|>
  filter(bdi_group != "NA")
#381

#want to see how many NAs there are for antidepressant column (6)
metadata_clean_small <- metadata_clean_small|>
  filter(antidepressant_on_off != "NA")
#375

##subsetting
#people with HIV, no HCV, depressed, use antidepressants
hiv_depressed_antidepressant <- metadata_clean_small|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == "On")
#39 people

#people with HIV, no HCV, depressed, no antidepressants
hiv_depressed_no_anti <- metadata_clean_small|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')|>
  filter(antidepressant_on_off == "Off")
#24 people in total

View(metadata_clean_edited)
#total number of people who has HIV and depression (and no HCV)
total_people <- metadata_clean_small|>
  filter(hiv_status_clean == 'HIV+')|>
  filter(hcv == 'NO')|>
  filter(bdi_group != 'minimal')
#63 people in total

#making the phyloseq
library(tidyverse)
library(phyloseq)
library(BiocManager)
library(microbiome)
library(indicspecies)

#Making the phyloseq
taxonomy = read.delim('taxonomy.tsv', row.names = 1)
counts = read.delim('feature-table.txt', skip=1, row.names=1)
tree = read_tree("tree.nwk")

otu_table <- otu_table(counts, taxa_are_rows = TRUE)

sample_data <- sample_data(total_people)

taxa_table <- taxonomy |>
  separate(col = Taxon, 
           into = c('Domain','Phylum','Class','Order',
                    'Family','Genus','Species'),
           sep=';', fill='right') |> 
  select(-Confidence)  |>
  as.matrix()

tax_table <- tax_table(taxa_table)

ps0 <- phyloseq(otu_table, sample_data, tax_table, tree)

saveRDS(ps0, "my_phyloseq_object_updated.rds")

ps = readRDS("my_phyloseq_object_updated.rds") |>
  tax_glom("Genus")

#ssubset phyloseq
on_anti = subset_samples(ps, antidepressant_on_off == 'On')
off_anti = subset_samples(ps, antidepressant_on_off == 'Off')

#core members
core_on_anti = core_members(on_anti, detection = 0.005, prevalence = 0.3)
core_off_anti = core_members(off_anti, detection = 0.005, prevalence = 0.3)

library(ggVennDiagram)

ggVennDiagram(list(core_on_anti, core_off_anti), set_size = 5, category.names = c('On', 'Off'))

#with filter
#indicator species analysis
ps_phylum = tax_glom(ps, 'Genus')
ps_relab = transform(ps_phylum, 'compositional')
ps_filt = filter_taxa(ps_relab, function(x) mean(x) > 0.001, TRUE) #keep that is true for
otu_table = data.frame(otu_table(ps_filt))

set.seed(400)
#indval = multipatt(t(otu_table), cluster = ps_filt@sam_data$antidepressant_on_off, control = how(nperm = 999))

indval = multipatt(t(otu_table), cluster = ps_filt@sam_data$antidepressant_on_off, control = how(nperm = 999))

summary(indval, indvalcomp = TRUE)
indval_table= as.data.frame(indval$sign)

View(indval_table)

#ind_val <- write.csv("indval_table")

#matching actual names to the codes
# extract the taxonomy table as a data frame
tax_info <- as.data.frame(tax_table(ps_filt))

# match the ISA results with the taxonomy names
indval_table$Genus <- tax_info[rownames(indval_table), "Genus"]

# move Genus to the first column so i can read
indval_table <- indval_table[, c("Genus", setdiff(names(indval_table), "Genus"))]

#plotting
phyla_to_plot = indval_table%>%
  filter(s.Off ==1 & s.On==0)%>%
  rownames()

df_of_taxa = prune_taxa(phyla_to_plot, ps_filt)%>%
  psmelt()

library(ggplot2)
plot_filter <- df_of_taxa%>%
  ggplot(aes(antidepressant_on_off, Abundance, fill=antidepressant_on_off))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=0, width=0.2)+
  facet_wrap(~Genus, ncol=4, scales = 'free')

#plot_filter + scale_fill_manual(values = c("Off" = "#6633FF","On"  = "#00FFCC" ))

library(scales)
final_plot <- plot_filter +
  scale_y_continuous(trans = pseudo_log_trans())+ scale_fill_manual(values = c("Off" = "#4876FF","On"  = "#FF82AB" )) 

final_plot

ggsave("final_plot.png", plot = final_plot, width = 10, height = 10, limitsize = FALSE)
