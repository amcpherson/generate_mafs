import dask.dataframe as dd


df = dd.read_csv(snakemake.input, sep='\t', comment='#', dtype='str')
df.to_csv(str(snakemake.output[0]), sep='\t', index=False, single_file=True)
