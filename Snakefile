import sys
import os
import time

from templates import *


# Config file
configfile: 'config/config.yaml'



########################################
#Collect information from Sample Sheet
########################################
run_check_bool, sample_data_dict, options_dict, sample_error_messages, options_error_messages, options_good_messages = messages(sample_sheet, fastq_dir)

# Sample IDs
sample_ids = sample_data_dict.keys()

# Make date directory
today = time.strftime('%Y%m%d')
date_dir = today + '/'
for d in config['directories']:
    os.makedirs(d, exist_ok=True)
os.makedirs(config['output_dir'] + date_dir, exist_ok=True)
os.makedirs(config['output_dir'] + date_dir + 'bam_files/', exist_ok=True)
os.makedirs(config['qc_temp'] + date_dir, exist_ok=True)

# Organism
organism_id = int(options_dict['organism'])
if organism_id==1:
    fa_link = config['mapping']['human']['download_link']['fasta']
    fa_file = config['mapping']['human']['fasta_file']
    gtf_link = config['mapping']['human']['download_link']['gtf']
    gtf_file = config['mapping']['human']['gtf_file']
elif organism_id==2:
    fa_link = config['mapping']['mouse']['download_link']['fasta']
    fa_file = config['mapping']['mouse']['fasta_file']
    gtf_link = config['mapping']['mouse']['download_link']['gtf']
    gtf_file = config['mapping']['mouse']['gtf_file']
elif organism_id==3:
    fa_link = options_dict['fa_link']
    gtf_link = options_dict['gtf_link']
    fa_file = options_dict['fa_link'].split('/')[-1]
    gtf_file = options_dict['gtf_link'].split('/')[-1]

# Extract genome ID
genome_id = fa_file.split('.')[1]

# Extract multimapping ID
multimap_n = int(options_dict['map_multimapper'])
if multimap_n == 1:
    multimap_id = 'noMM'
else:
    multimap_id = 'MM'

# Extract counting parameter
strand_id = int(options_dict['count_strandedness'])
if strand_id == 1:
    strand_mode = 'reverse'
elif strand_id == 2:
    strand_mode = 'yes'
elif strand_id == 3:
    strand_mode = 'no'
align_id = int(options_dict['count_align_mode'])
if align_id == 1:
    align_mode = 'intersection-strict'
elif align_id == 2:
    align_mode = 'union'
elif align_id == 3:
    align_mode = 'intersection-nonempty'

# Extract bigwig bin size
bw_bin_size = int(options_dict['bigwig'])


#######################################
# Workflow
#######################################

# Run only if Sample Sheet filled correct
if run_check_bool:
    print(options_good_messages)
    # Proces fastq files from input folder (input/fastq/) or SRA Archive
    if options_dict['seq_type'] == '1': # Paired end reads
        if options_dict['input_type'] == '1':  # From uploaded FASTQ files
            include: 'rules/01_rename_fq.smk'
        if options_dict['input_type'] == '2':  # From SRA Archive
            include: 'rules/01_fastq_download_paired_end.smk'
        include: 'rules/02_fastq_quality_control_paired.smk'
        if options_dict['trim_detection'] == '1': # Adapters auto-detected
            include: 'rules/03_trim_fastq_paired_autodetect_adapters.smk'
        else: # Special adapters manually trimmed
            include: 'rules/03_trim_fastq_paired_manual_adapters.smk'
        include: 'rules/04_extract_genome_files.smk'
        include: 'rules/04_index_build_star.smk'
        include: 'rules/05_align_star_paired.smk'
        include: 'rules/06_count_htseq_paired.smk'
        if len(sample_ids) > 1:
            include: 'rules/06_merge_counts.smk'
            include: 'rules/07_normalise_counts.smk'
        include: 'rules/08_make_bigwigs_paired.smk'
        include: 'rules/09_multiqc_report_paired.smk'
        include: 'rules/10_save_pipeline_parameters.smk'

    elif options_dict['seq_type'] == '2': # Single reads
        if options_dict['input_type'] == '1':  # From uploaded FASTQ files
            include: 'rules/01_rename_fq.smk'
        if options_dict['input_type'] == '2':  # From SRA Archive
            include: 'rules/01_fastq_download_single_end.smk'
        include: 'rules/02_fastq_quality_control_single.smk'
        if options_dict['trim_detection'] == '1': # Adapters auto-detected
            include: 'rules/03_trim_fastq_single_autodetect_adapters.smk'
        else:  # Special adapters manually trimmed -
            include: 'rules/03_trim_fastq_single_manual_adapters.smk'
        include: 'rules/04_extract_genome_files.smk'
        include: 'rules/04_index_build_star.smk'
        include: 'rules/05_align_star_single.smk'
        include: 'rules/06_count_htseq_single.smk'
        if len(sample_ids) > 1:
            include: 'rules/06_merge_counts.smk'
            include: 'rules/07_normalise_counts.smk'
        include: 'rules/08_make_bigwigs_single.smk'
        include: 'rules/09_multiqc_report_single.smk'
        include: 'rules/10_save_pipeline_parameters.smk'
else:
    print(options_error_messages)
    print(sample_error_messages)


#######################################
# Execute pipeline
#######################################

# Run only if Sample Sheet filled correct
if run_check_bool:
    if len(sample_ids) > 1:
        if options_dict['seq_type'] == '1':
            rule all:
                input:
                    expand(config['output_dir'] + date_dir + 'bam_files/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam.bai', sample_id=sample_ids),
                    config['output_dir'] + date_dir + 'counts/counts_normalised_deseq2.csv',
                    config['output_dir'] + date_dir + 'counts/counts_normalised_tpm.csv',
                    expand(config['output_dir'] + date_dir + 'bigwig/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP_{norm_id}_bin' + str(bw_bin_size) +'_{strand}.bw',
                        sample_id=sample_ids, norm_id=['noNorm','sfNorm'], strand=['plus','minus']),
                    config['output_dir'] + date_dir + 'multiqc_report.html',
                    config['output_dir'] + date_dir + 'pipeline_parameters/templates.py'
        else:
            rule all:
                input:
                    expand(config['output_dir'] + date_dir + 'bam_files/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam.bai',sample_id=sample_ids),
                    config['output_dir'] + date_dir + 'counts/counts_normalised_deseq2.csv',
                    config['output_dir'] + date_dir + 'counts/counts_normalised_tpm.csv',
                    expand(config[ 'output_dir'] + date_dir + 'bigwig/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP_{norm_id}_bin' + str(bw_bin_size) + '.bw',
                        sample_id=sample_ids,norm_id=['noNorm', 'sfNorm']),
                    config['output_dir'] + date_dir + 'multiqc_report.html',
                    config['output_dir'] + date_dir + 'pipeline_parameters/templates.py'

    else:
        if options_dict['seq_type'] == '1':
            rule all:
                input:
                    expand(config['output_dir'] + date_dir + 'bam_files/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam.bai',sample_id=sample_ids),
                    expand(config['output_dir'] + date_dir + 'bigwig/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP_noNorm_bin' + str(bw_bin_size) +'_{strand}.bw',
                       sample_id=sample_ids, strand=['plus', 'minus']),
                    config['output_dir'] + date_dir + 'multiqc_report.html',
                    config['output_dir'] + date_dir + 'pipeline_parameters/templates.py'
        else:
            rule all:
                input:
                    expand(config['output_dir'] + date_dir + 'bam_files/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam.bai',sample_id=sample_ids),
                    expand(config['output_dir'] + date_dir + 'bigwig/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP_noNorm_bin' + str(bw_bin_size) + '.bw',
                        sample_id=sample_ids),
                    config['output_dir'] + date_dir + 'multiqc_report.html',
                    config['output_dir'] + date_dir + 'pipeline_parameters/templates.py'



