library(tidyverse)
library(ANCOMBC)
library(phyloseq)
library(ggplot2)

ps = readRDS("my_phyloseq_object.rds") |>
  tax_glom("Genus")

