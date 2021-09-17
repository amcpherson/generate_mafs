import dask.dataframe as dd


df = dd.read_csv(snakemake.input[0], sep='\t', comment='#', dtype='str')
df = df[df['ONCOGENIC'].isin(['Likely Oncogenic', 'Predicted Oncogenic', 'Oncogenic'])]
df.to_csv(str(snakemake.output[0]), sep='\t', index=False, single_file=True)

