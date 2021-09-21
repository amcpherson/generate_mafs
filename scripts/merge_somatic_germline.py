import dask.dataframe as dd

print(dir(snakemake))

somatic = dd.read_csv(snakemake.input[0], sep='\t', comment='#', dtype='str').assign(Mutation_Status="Somatic")
germline = dd.read_csv(snakemake.input[1], sep='\t', comment='#', dtype='str').assign(Mutation_Status="Germline")
dd.concat([somatic, germline]).to_csv(str(snakemake.output[0]), sep='\t', index=False, single_file=True)
