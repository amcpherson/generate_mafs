library(data.table)
library(tidyverse)
library(vroom)

col_types <- cols(
  t_depth = col_integer(),
  t_alt_count = col_integer(),
  t_ref_count = col_integer(),
  n_depth = col_integer(),
  n_alt_count = col_integer(),
  n_ref_count = col_integer(),
  Start_Position = col_integer(),
  End_Position = col_integer(),
  gnomAD_AF = col_double(),
  .default = col_character())
maf <- vroom(snakemake@input[["somatic"]], comment="#", col_types=col_types)
maf_germline <- vroom(snakemake@input[["germline"]], comment="#", col_types=col_types)

#filter based on recurrent mutations across tumors
maf_filt <- maf %>%
  select(-FILTER) %>%
  group_by(across(c(-t_depth, -t_alt_count,-t_ref_count,-n_depth, -n_alt_count,-n_ref_count))) %>%
  # sum counts across libraries
  summarise(t_depth = sum(t_depth, na.rm = T), 
            t_alt_count = sum(t_alt_count, na.rm = T),
            t_ref_count = sum(t_ref_count, na.rm = T),
            n_depth = sum(n_depth, na.rm = T), 
            n_alt_count = sum(n_alt_count, na.rm = T),
            n_ref_count = sum(n_ref_count, na.rm = T)) %>%
  ungroup() %>%
  #count number of times the same variant appears in different tumors
  group_by(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(nsample = length(unique(Tumor_Sample_Barcode))) %>% 
  ungroup() %>% 
  # remove variants that appear > 1 tumor but retain TP53 muts
#  filter(!(nsample > 1 & 
#             SYMBOL != "TP53")) %>%
  select(-nsample) %>%
  mutate(VAF = t_alt_count / t_depth) %>%
  mutate(Mutation_Status = "Somatic")

#filter based on recurrent mutations across tumors
maf_germline_filt <- maf_germline %>%
  select(-FILTER) %>%
  group_by(across(c(-t_depth, -t_alt_count,-t_ref_count,-n_depth, -n_alt_count,-n_ref_count))) %>%
  # sum counts across libraries
  summarise(t_depth = sum(t_depth, na.rm = T), 
            t_alt_count = sum(t_alt_count, na.rm = T),
            t_ref_count = sum(t_ref_count, na.rm = T),
            n_depth = sum(n_depth, na.rm = T), 
            n_alt_count = sum(n_alt_count, na.rm = T),
            n_ref_count = sum(n_ref_count, na.rm = T)) %>%
  ungroup() %>%
  #count number of times the same variant appears in different tumors
  group_by(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(nsample = length(unique(Tumor_Sample_Barcode))) %>% 
  ungroup() %>% 
  # remove variants that appear > 1 tumor but retain TP53 muts
#  filter(!(nsample > 1 & 
#             !(SYMBOL %in% c("TP53", "BRCA1", "BRCA2")))) %>%
  # filter for oncogenic
  filter(ONCOGENIC %in% c("Oncogenic", "Likely Oncogenic")) %>%
  select(-nsample) %>%
  mutate(VAF = t_alt_count / t_depth) %>%
  mutate(Mutation_Status = "Germline")

#filter based on population AF, mutation impact/consequence
maf_filt <- maf_filt %>% 
  #filter for gnomad AF of 1%
  filter(gnomAD_AF < 0.01 | is.na(gnomAD_AF)) %>% 
  #filter for protein coding or Oncogenic
  filter(!Variant_Classification %in% c("Silent", "Intron", "IGR", "3'UTR", "5'UTR", "3'Flank", "5'Flank")) %>%
  #filter for high impact or Oncogenic
  filter(IMPACT %in% c("HIGH", "MODERATE", "LOW") | ONCOGENIC != "")

#filter based on population AF, mutation impact/consequence
maf_germline_filt <- maf_germline_filt %>% 
  #filter for gnomad AF of 1%
  filter(gnomAD_AF < 0.01 | is.na(gnomAD_AF))

#merge germline and somatic (remove vcf id column due to type mismatch)
maf_merged <- bind_rows(maf_filt %>% select(-vcf_id, -vcf_qual), maf_germline_filt %>% select(-vcf_id, -vcf_qual))

fwrite(maf_merged, file = snakemake@output[[1]], sep = "\t")
  
