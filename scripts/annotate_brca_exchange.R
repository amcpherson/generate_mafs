library(data.table)
library(tidyverse)

brcaex <- fread(snakemake@config$brcaexchange) %>% 
  separate(genomic_vcf37_source, c("chr_hg19", "pos_hg19", "ref_alt_hg19"), remove = FALSE, sep = ":") %>% 
  mutate(chr_hg19 = str_remove(chr_hg19, "chr"),
         pos_hg19 = as.numeric(str_remove(pos_hg19, "g."))) %>%
  separate(ref_alt_hg19, c("ref_hg19", "alt_hg19"), sep = ">") %>% 
  select(genomic_vcf37_source, chr_hg19, pos_hg19, ref_hg19, alt_hg19, clinical_significance_enigma) %>% 
  mutate(brcaexchange = TRUE) %>% 
  rename(clinical_significance_brcaexchange = clinical_significance_enigma)
           
maf <- fread(snakemake@input[[1]])

merged <- brcaex[maf, on = .(chr_hg19 == Chromosome, pos_hg19 == Start_Position, ref_hg19 == Reference_Allele, alt_hg19 == Tumor_Seq_Allele2)]

fwrite(merged, snakemake@output[[1]], sep = "\t")
