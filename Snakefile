# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

configfile: 'config.yaml'

import os
import yaml
import pandas as pd

metadata = pd.read_csv(config['metadata'])
oncokb_api_key = config['oncokb_api_key']
intermediate_dir = config['intermediate_dir']
results_dir = config['results_dir']
scripts_dir = config['scripts_dir']

# # List of unique tumour samples for somatic calls
# tumor_sample_ids = metadata.query('target_sample_category != "NORMAL"')['target_sample'].unique()
#
# # List of unique normal samples for germline calls
# normal_sample_ids = metadata.query('target_sample_category == "NORMAL"')['target_sample'].unique()
#
# # List of paths to input files
# paths = metadata.set_index(['target_sample', 'app', 'file_type'])['path'].to_dict()
#
# # We need a set of germline calls for every tumour sample, thus
# # we need to duplicate germline calls where there are multiple
# # tumour samples per patient.  Generate a list of tumour sample ids
# # for each normal sample id to allow replication of germline calls
# # for each tumour sample.
# normal_tumor_map = {}
# for individual, indiv_metadata in metadata.groupby('individual'):
#     indiv_normal = indiv_metadata.query('target_sample_category == "NORMAL"')['target_sample'].unique()
#     indiv_tumor = indiv_metadata.query('target_sample_category != "NORMAL"')['target_sample'].unique()
#     for normal_sample_id in indiv_normal:
#         for tumor_sample_id in indiv_tumor:
#             if normal_sample_id not in normal_tumor_map:
#                 normal_tumor_map[normal_sample_id] = []
#             normal_tumor_map[normal_sample_id].append(tumor_sample_id)
#
# # Generate a mapping of normal sample ids used for each tumour sample
# # during somatic calling.  We will add these sample ids into the maf,
# # replacing existing sample ids
# tumour_normal_map = metadata.query('target_sample_category != "NORMAL"').set_index(['target_sample'])['reference_sample'].to_dict()
#
rule all:
    input:
        os.path.join(results_dir, 'cohort.maf'),
        os.path.join(results_dir, 'cohort_filtered.maf')

# Somatic
#

# def get_somatic_input_paths(wildcards):
#     return {
#         'consensus_somatic_maf': paths[(wildcards.tumor_sample_id, 'WGS-SOMATICCALLING', 'consensus_somatic_maf')],
#         'museq_paired_annotated': paths[(wildcards.tumor_sample_id, 'WGS-SOMATICCALLING', 'museq_paired_annotated')],
#         'strelka_snv_annotated': paths[(wildcards.tumor_sample_id, 'WGS-SOMATICCALLING', 'strelka_snv_annotated')],
#     }
#
# rule somatic_filter_maf:
#     input: unpack(get_somatic_input_paths)
#     output: os.path.join(intermediate_dir, 'somatic_{tumor_sample_id}.filtered.maf')
#     params: scripts_dir=scripts_dir
#     singularity: "docker://amcpherson/filtermafs"
#     shell: 'python {params.scripts_dir}/filter_snv_vcf.py {input.consensus_somatic_maf} {input.museq_paired_annotated} {input.strelka_snv_annotated} {output}'
#
# rule somatic_fix_sample_id:
#     input: os.path.join(intermediate_dir, 'somatic_{tumor_sample_id}.filtered.maf')
#     output: os.path.join(intermediate_dir, 'somatic_{tumor_sample_id}.filtered.fixup.maf')
#     params:
#         tumor_sample_ids = lambda wildcards: [wildcards.tumor_sample_id],
#         normal_sample_id = lambda wildcards: tumour_normal_map[wildcards.tumor_sample_id],
#     singularity: "docker://amcpherson/filtermafs"
#     script: "scripts/fix_sample_id.py"
#
# rule somatic_annotate_maf:
#     input: os.path.join(intermediate_dir, 'somatic_{tumor_sample_id}.filtered.fixup.maf')
#     output: os.path.join(intermediate_dir, 'somatic_{tumor_sample_id}.filtered.fixup.annotated.maf')
#     singularity: "docker://amcpherson/oncokb-annotator"
#     shell: 'python /oncokb-annotator/MafAnnotator.py -i {input} -o {output} -b {oncokb_api_key}'
#
# rule somatic_merge_mafs:
#     input: expand(os.path.join(intermediate_dir, 'somatic_{tumor_sample_id}.filtered.fixup.annotated.maf'), tumor_sample_id=tumor_sample_ids)
#     output: os.path.join(intermediate_dir, 'somatic.maf')
#     singularity: "docker://amcpherson/filtermafs"
#     script: "scripts/merge_mafs.py"

# Germline
# 

# def get_germline_input_paths(wildcards):
#     return [paths[(wildcards.normal_sample_id, 'WGS-GERMLINECALLING', 'consensus_germline_maf')]]
#
# rule germline_filter_maf:
#     input: unpack(get_germline_input_paths)
#     output: os.path.join(intermediate_dir, 'germline_{normal_sample_id}.filtered.maf')
#     singularity: "docker://amcpherson/filtermafs"
#     script: "scripts/filter_germline_maf.py"
#
# rule germline_fix_sample_id:
#     input: os.path.join(intermediate_dir, 'germline_{normal_sample_id}.filtered.maf')
#     output: os.path.join(intermediate_dir, 'germline_{normal_sample_id}.filtered.fixup.maf')
#     params:
#         normal_sample_id = lambda wildcards: wildcards.normal_sample_id,
#         tumor_sample_ids = lambda wildcards: normal_tumor_map.get(wildcards.normal_sample_id, []),
#     singularity: "docker://amcpherson/filtermafs"
#     script: "scripts/fix_sample_id.py"
#
# rule germline_annotate_maf:
#     input: os.path.join(intermediate_dir, 'germline_{normal_sample_id}.filtered.fixup.maf')
#     output: os.path.join(intermediate_dir, 'germline_{normal_sample_id}.filtered.fixup.annotated.maf')
#     singularity: "docker://amcpherson/oncokb-annotator"
#     shell: 'python /oncokb-annotator/MafAnnotator.py -i {input} -o {output} -b {oncokb_api_key}'
#
# rule germline_merge_mafs:
#     input: expand(os.path.join(intermediate_dir, 'germline_{normal_sample_id}.filtered.fixup.annotated.maf'), normal_sample_id=normal_sample_ids)
#     output: os.path.join(intermediate_dir, 'germline.maf')
#     singularity: "docker://amcpherson/filtermafs"
#     script: "scripts/merge_mafs.py"
#
# Cohort somatic and germline
#

# rule merge_somatic_germline:
#     input:
#         os.path.join(intermediate_dir, 'somatic.maf'),
#         os.path.join(intermediate_dir, 'germline.maf')
#     output: os.path.join(intermediate_dir, 'merged.maf')
#     singularity: "docker://amcpherson/filtermafs"
#     script: "scripts/merge_somatic_germline.py"

rule annotate_brcaexchange:
    #input: os.path.join(intermediate_dir, 'merged.maf')
    input: os.path.join(intermediate_dir, 'germline.maf')
    output: os.path.join(results_dir, 'cohort.maf')
    singularity: "docker://amcpherson/filtermafs"
    script: "scripts/annotate_brca_exchange.py"

rule filter_maf:
    input: os.path.join(results_dir, 'cohort.maf')
    output: os.path.join(results_dir, 'cohort_filtered.maf')
    singularity: "docker://amcpherson/filtermafs"
    script: "scripts/filter_maf.py"

