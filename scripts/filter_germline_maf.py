import dask.dataframe as dd


df = dd.read_csv(snakemake.input, sep='\t', comment='#', dtype='str')

exclude_classes = [
    "Silent",
    "Intron",
    "IGR",
    "3'UTR",
    "5'UTR",
    "3'Flank",
    "5'Flank",
]

df = df[~df['Variant_Classification'].isin(exclude_classes)]

df['gnomAD_AF'] = df['gnomAD_AF'].astype(float)
df = df.query('gnomAD_AF.isnull() or gnomAD_AF < 0.01')

df['SIFT_call'] = df['SIFT'].str.replace('\(.*', '', regex=True)
df = df.query('SIFT_call != "tolerated"')

df['PolyPhen_call'] = df['PolyPhen'].str.replace('\(.*', '', regex=True)
df = df.query('PolyPhen_call != "benign"')

df.to_csv(str(snakemake.output), sep='\t', index=False, single_file=True)
