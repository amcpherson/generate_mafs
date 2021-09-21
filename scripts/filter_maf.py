import dask.dataframe as dd


df = dd.read_csv(snakemake.input, sep='\t', comment='#', dtype='str')

# Filter all non-oncogenic with the exception of TP53
oncogenic = ['Likely Oncogenic', 'Predicted Oncogenic', 'Oncogenic', 'Inconclusive']
df = df.query(f'Hugo_Symbol == "TP53" | ONCOGENIC in {oncogenic}')

# Filter all BRCA that are Benign (note only BRCA will have a valid value in this column)
df = df.query('clinical_significance_brcaexchange != "Benign"')

df.to_csv(str(snakemake.output), sep='\t', index=False, single_file=True)

