import dask.dataframe as dd


df = dd.read_csv(snakemake.input, sep='\t', comment='#', dtype='str')

fixed_df = []
for tumor_sample_id in snakemake.params['tumor_sample_ids']:
    df['Tumor_Sample_Barcode'] = tumor_sample_id
    df['Matched_Norm_Sample_Barcode'] = snakemake.params['normal_sample_id']
    fixed_df.append(df.copy())
fixed_df = dd.concat(fixed_df, ignore_index=True)

fixed_df.to_csv(str(snakemake.output[0]), sep='\t', index=False, single_file=True)
