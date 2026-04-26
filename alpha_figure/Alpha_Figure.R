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

# load data 
metadata = read.csv("metadata_filtered.csv")
metadata$antidepressant_on_off = as.factor(metadata$antidepressant_on_off)

metadata = metadata %>% 
  mutate(antidepressant_on_off = recode(antidepressant_on_off, 
                                        "0" = "Off",
                                        "1" = "On")) %>% 
  filter(!is.na(antidepressant_on_off))

hiv_depressed_antidepressant <- metadata  %>% 
  filter(antidepressant_on_off == 'On') # 38

hiv_depressed_no_anti <- metadata %>% 
  filter(antidepressant_on_off == 'Off') # 23

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_id

metadata <- metadata %>% 
  select(-X, -KEEP., -WHY., -host_subject_id, -Collection_Date...4, -antidepressant_count...5)

write.table(metadata, file = "metadata_filtered.tsv", sep = "\t", row.names = FALSE)

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
ps <- phyloseq(sample_data(metadata),
               otu_table(counts_formatted, taxa_are_rows = T),
               tax_table(taxonomy_formatted))

hist(sample_sums(ps))
rarecurve(t(data.frame(ps@otu_table)),
          step=1000, 
          label = FALSE)

psrare = ps %>% 
  rarefy_even_depth(sample.size = 10000, rngseed = 316)
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
  scale_fill_manual(values = c("Off" = "#4876FF","On"  = "#FF82AB")) + 
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
