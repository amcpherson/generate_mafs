library(data.table)
library(tidyverse)
library(vroom)

maf <- vroom(snakemake@input[[1]], comment="#")

#filter based on population AF, mutation impact/consequence
maf_filt <- maf %>%
  #filter common_variant in 
  filter(grepl('common_variant', FILTER)) %>%
  #filter for AF of 1%
  filter(AF < 0.01 | is.na(AF)) %>% 
  #filter for gnomad AF of 1%
  filter(gnomAD_AF < 0.01 | is.na(gnomAD_AF)) %>% 
  #filter for protein coding or Oncogenic
  filter(!Variant_Classification %in% c("Silent", "Intron", "IGR", "3'UTR", "5'UTR", "3'Flank", "5'Flank")) %>%
  #filter for impact
  filter(IMPACT %in% c("HIGH", "MODERATE", "LOW"))

fwrite(maf_filt, file = snakemake@output[[1]], sep = "\t")
