# load packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
library(microbiome)
library(indicspecies)
library(writexl)
library(ANCOMBC)
library(readxl)

# loading and tidying dataset
metadata = read.delim('depression_metadata_full.txt', row.names=1)  

metadata$antidepressant_on_off = as.factor(metadata$antidepressant_on_off)

metadata_clean <- metadata %>%
  select(antidepressant_on_off, antidepressant_count, bdi_group, bdi_ii, hiv_status_clean, hcv) %>%
  drop_na(antidepressant_on_off, antidepressant_count, bdi_group, bdi_ii, hiv_status_clean, hcv) %>%
  filter(hcv == "NO")|>
  filter(hiv_status_clean == "HIV+") %>%
  filter(bdi_group != 'minimal')

data = metadata_clean %>% 
  mutate(antidepressant_on_off = recode(antidepressant_on_off, 
                                        "0" = "Off",
                                        "1" = "On")) %>% 
  filter(!is.na(antidepressant_on_off))

hiv_depressed_antidepressant <- data  %>% 
  filter(antidepressant_on_off == 'On')
#78 people in total

hiv_depressed_no_anti <- data %>% 
  filter(antidepressant_on_off == 'Off')
#51 people in total

# loading and tidying taxonomy and counts
# Taxonomy #
taxonomy = read.delim('taxonomy.tsv', row.names = 1) 

taxonomy_formatted = taxonomy %>% 
  separate(col = Taxon,
           into = c("Domain", "Phylum","Class", "Order",
                    "Family", "Genus", "Species"),
           sep = ";", fill = "right") %>% 
  select(-Confidence) %>% 
  as.matrix()
# Count #
counts = read.delim('feature-table.txt', skip=1, row.names=1) 

counts_formatted = counts %>% 
  as.matrix()

# Figure Generation
ps <- phyloseq(sample_data(data),
               otu_table(counts_formatted, taxa_are_rows = T),
               tax_table(taxonomy_formatted))

hist(sample_sums(ps))
rarecurve(t(data.frame(ps@otu_table)),
          step=1000, 
          label = FALSE)

psrare = ps %>% 
  rarefy_even_depth(sample.size = 20000, rngseed = 316)
table(sample_sums(psrare))

set.seed(315)

p = plot_richness(psrare, x = "antidepressant_on_off",
                  measures = c("Observed", "Shannon", "Chao1", "Simpson"),
                  color = "antidepressant_on_off")
p

pdata = p$data
str(pdata)

p_formatted = pdata %>% 
  ggplot(aes(antidepressant_on_off,value,
             fill = antidepressant_on_off)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height=0,width=0.2,size = 0.5, alpha = 0.5) +
  theme_classic(base_size=14) +
  facet_wrap(~variable,ncol=4,scales = 'free_y') +
  ylab('Alpha Diversity') + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1)) +
  labs(fill='Antidepressant Use') + 
  scale_fill_grey(start=0.6, end=1) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_formatted

set.seed(316) 
comparisons = list(c("Off", "On"))
final_figure = p_formatted + 
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     size = 4) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
final_figure
ggsave("alpha_diversity.png", plot = final_figure)
