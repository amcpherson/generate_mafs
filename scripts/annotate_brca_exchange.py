import dask.dataframe as dd
import pandas as pd
import dask
dask.config.set(scheduler='synchronous')


brca_exchange = pd.read_csv(snakemake.config['brcaexchange'], sep='\t')

brca_exchange['Chromosome'] = brca_exchange['genomic_vcf37_source'].str.split(':', expand=True)[0].str.replace('chr', '', regex=False)
brca_exchange['Start_Position'] = brca_exchange['genomic_vcf37_source'].str.split(':', expand=True)[1].str.replace('g.', '', regex=False).astype(int)
brca_exchange['ref_alt'] = brca_exchange['genomic_vcf37_source'].str.split(':', expand=True)[2]
brca_exchange['Reference_Allele'] = brca_exchange['ref_alt'].str.split('>', expand=True)[0]
brca_exchange['Tumor_Seq_Allele2'] = brca_exchange['ref_alt'].str.split('>', expand=True)[1]
brca_exchange['clinical_significance_brcaexchange'] = brca_exchange['clinical_significance_enigma']

brca_exchange = brca_exchange[[
    'Chromosome',
    'Start_Position',
    'Reference_Allele',
    'Tumor_Seq_Allele2',
    'clinical_significance_brcaexchange',
]]

df = pd.read_csv(snakemake.input[0], sep='\t', comment='#', dtype='str')

df['Start_Position'] = df['Start_Position'].astype(int)

df = df.merge(brca_exchange, how='left')

df.to_csv(str(snakemake.output), sep='\t', index=False, single_file=True)

