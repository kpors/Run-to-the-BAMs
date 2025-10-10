# Run-To-The-BAMs

##  - An automatic and efficient Snakemake workflow for RNA-seq data processing - 

[TO BE UPDATED, KPJ 20251010]

## Pipeline overview

![](https://github.com/kpors/Run-to-the-BAMs/misc/rna_seq_snakemake.jpeg)

The Run-To-The_BAMs pipeline automates RNA-seq data processing from raw FASTQ files to quality-checked, mapped, and normalised gene counts, along with visualisation-ready BigWig files.

In this pipeline, you either upload our own FASTQ files or select SRR IDs from the SRA Archive for analysis. The pipeline has following steps:

- Quality control (FastQC, Trim Galore!, MultiQC)
  - Output: Read and mapping quality control report
- Mapping reads to selected genome (STAR)
  - Output: BAM and BAI file
- Count reads mapped to genome (HTSeq)
  - Output: Raw gene counts in CSV format including gene id, gene name and gene biotype
- Normalise counts by DeSeq2 and TPM methods (DeSeq2, R)
  - Output: Normalised gene counts in CSV format including gene id, gene name and gene biotype
- Generate BigWig file (deepTools)
  - Output: Normalised and raw BigWig files for both plus and minus strand



## Requirements

This pipeline is designed to run on a **SLURM-based high-performance computing (HPC) cluster**.  
It uses **Snakemake** to submit and manage jobs automatically through SLURM.  

To efficiently process RNA-seq data, the HPC node should have at least **64 GB of RAM** and **16 CPU cores**, but when running several samples or workflow rules in parallel, proportionally higher resources are required.



## Installation

#### Run-to-the-bams

On your HPC cluster, navigate to the directory designated for storing the pipeline files before downloading them. 

Download pipeline files:

```
bash

$ git clone https://github.com/kpors/Run-to-the-BAMs
```



#### Conda

In this manual, we use **Conda** to manage the software environment required for running the pipeline.  
It is therefore important that you have Conda installed.  
If Conda is not already available on your system, you can install it by following the commands below, line by line.

##### Conda installation

```
$ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh
$ chmod +x miniforge.sh
$ bash miniforge.sh -b
$ ./miniforge3/bin/conda init bash
```

##### Conda configuration

These settings ensure that Conda installs packages from the correct bioinformatics channels with strict version control.

```
$ conda config --append channels bioconda
$ conda config --append channels genomedk
$ conda config --set channel_priority strict
$ conda config --set auto_activate_base false
```



#### Snakemake installation

To enable the use of Snakemake and its associated tools, install it in a Conda environment.

```
$ conda create -n snakemake snakemake-executor-plugin-slurm snakemake
```

Activate the environment 

```
$ conda activate snakemake
```



#### Snakemake pipeline execution

##### Processing own FASTQ files

Upload FASTQ files to Run-to-the-bams/input/fastq/

Fill out the sample sheet in Run-to-the-bams/input/ and save it as a CSV exactly as this:

Run-to-the-bams/input/Sample_sheet.csv

##### Processing FASTQ files from SRA Archive (https://www.ncbi.nlm.nih.gov/sra)

Fill out the sample sheet in Run-to-the-bams/input/ with SRRXXXXXXX IDs and save it as a CSV exactly as this:

Run-to-the-bams/input/Sample_sheet.csv

Run the pipeline

Test if everything is alright (dry-run):

```
(snakemake) $ snakemake --wait-for-files --executor slurm --workflow-profile profiles/slurm --use-conda --jobs unlimited --dry-run 
```

Run to the BAMs:

```
(snakemake) $ snakemake --wait-for-files --executor slurm --workflow-profile profiles/slurm --use-conda --jobs unlimited
```

After the pipeline is done output files can be found in the output directory:
Run-to-the-bams/output/
