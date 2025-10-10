# Run To The BAMs
## - A Snakemake pipeline to process RNA sequencing data

SAMPLE ID & FILE NAMES
<sample_id>_R1.fastq.gz
<sample_id>_R2.fastq.gz
<SRA_id>_R1.fastq.gz
<SRA_id>_R2.fastq.gz

DATA INPUT
For each run it is only possible to process fastq files from
    - either input folder (input/fastq/) or SRA archive
AND only possible to process fastq data type that are
    - either paired-end or single read
AND only reads with the same read length

If you need to proces both fastq files from input folder and SRA archive or
fastq data that are sequenced as paired-end and single reads or with different read lengths, please run these separately



SAMPLE SHEET
It is RECOMMENDED to open and edit the sample sheet in either Excel or Numbers.
It is OBLIGATORY to save/export the sample sheet as CSV and put it in the input folder


TRIMMING
By default adapters are auto-detected with Trim-Galore, that is used to trim the fastq files. If you want to put in
adapter sequences manually, you need to do following changes in the ./config/config.yaml file:
    - Change ['trim_fastq']['special_adapters']['bool'] to False
    - Add your first adapter sequence in ['trim_fastq']['special_adapters']['adapter_1']
    - Add your second adapter sequence in ['trim_fastq']['special_adapters']['adapter_2']

