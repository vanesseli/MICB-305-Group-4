#phyloseq object creation

library(tidyverse)
library(phyloseq)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq", force = TRUE)

#Loading in the stuff


#Making the phyloseq
