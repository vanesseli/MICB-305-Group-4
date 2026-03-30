library(tidyverse)
library(phyloseq)
library(microbiome)
library(indicspecies)

set.seed(400)

ps = readRDS("../project_analyses/my_phyloseq_object.rds") %>%
  tax_glom(ps,'Genus')
ps_relab = transform(ps, 'compositional')
ps_filt = filter_taxa(ps_relab, function(x) mean(x) > 0.001, TRUE)
otu_table = data.frame(otu_table(ps_filt))

#indval = multipatt(t(otu_table), cluster = ps_filt@sam_data$antidepressant_on_off, control = how(nperm = 999))

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

indval_table = indval_table %>%
  filter(p.value < 0.05)

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

final_plot <- plot_filter +
  scale_y_log10()+ scale_fill_manual(values = c("Off" = "#6633FF","On"  = "#00FFCC" ))

final_plot
