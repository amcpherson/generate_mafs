import pandas as pd
import vcf
import fire


def read_vcf(vcf_filename):
    """ Read a vcf and return basic SNV data.
    """

    vcf_reader = vcf.Reader(open(vcf_filename, 'rb'))

    data = []
    for record in vcf_reader:
        for alt in record.ALT:
            data.append((str(record.CHROM), int(record.POS), str(record.REF), str(alt), record.QUAL))
    data = pd.DataFrame(data, columns=['chrom', 'coord', 'ref', 'alt', 'qual'])

    return data


exclude_classes = [
    "Silent",
    "Intron",
    "IGR",
    "3'UTR",
    "5'UTR",
    "3'Flank",
    "5'Flank",
]


def filter_vcf(consensus_somatic_maf, museq_paired_annotated, strelka_snv_annotated, output_maf):
    """ Filter maf using museq and strelka vcf.
    """

    wgs_maf_data = pd.read_csv(
        consensus_somatic_maf,
        sep='\t', comment='#', dtype=str, low_memory=False)

    wgs_maf_data['chrom'] = wgs_maf_data['Chromosome']
    wgs_maf_data['coord'] = wgs_maf_data['Start_Position'].astype(int)
    wgs_maf_data['ref'] = wgs_maf_data['Reference_Allele']
    wgs_maf_data['alt'] = wgs_maf_data['Tumor_Seq_Allele2']

    museq_vcf = read_vcf(museq_paired_annotated)
    strelka_vcf = read_vcf(strelka_snv_annotated)

    wgs_maf_data = wgs_maf_data.merge(
        museq_vcf.assign(museq_snv_predicted=1),
        how='left', on=['chrom', 'coord', 'ref', 'alt'],
        suffixes=('', '_museq_snv'))

    wgs_maf_data = wgs_maf_data.merge(
        strelka_vcf.assign(strelka_snv_predicted=1),
        how='left', on=['chrom', 'coord', 'ref', 'alt'],
        suffixes=('', '_strelka_snv'))

    wgs_maf_data['museq_snv_predicted'] = wgs_maf_data['museq_snv_predicted'].fillna(0).astype(int)
    wgs_maf_data['strelka_snv_predicted'] = wgs_maf_data['strelka_snv_predicted'].fillna(0).astype(int)

    # Filter SNPs not predicted by both callers
    wgs_maf_data = wgs_maf_data.query('Variant_Type != "SNP" or (museq_snv_predicted == 1 and strelka_snv_predicted == 1)')

    wgs_maf_data = wgs_maf_data.drop(['chrom', 'coord', 'ref', 'alt', 'qual', 'museq_snv_predicted', 'qual_strelka_snv', 'strelka_snv_predicted'], axis=1)

    wgs_maf_data = wgs_maf_data[~wgs_maf_data['Variant_Classification'].isin(exclude_classes)]

    wgs_maf_data.to_csv(output_maf, sep='\t', index=False)


if __name__ == '__main__':
    fire.Fire(filter_vcf)


