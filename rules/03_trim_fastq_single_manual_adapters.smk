rule trim_fastq_single_manual:
    # Description
    input:
        fq1=lambda wc: config['fastq_dir'] + sample_data_dict[wc.sample_id][0].replace('fq.gz','fastq.gz')
    output:
        config['temp_fastq'] + '{sample_id}_trimmed.fq'
    conda:
        '../' + config['envs_dir'] + 'trim_galore.yml'
    threads: config['slurm']['trim_fastq']['threads']
    resources:
        mem=config['slurm']['trim_fastq']['mem'],
        runtime=config['slurm']['trim_fastq']['runtime']
    params:
        outdir_fastqc = config['qc_temp'] + date_dir,
        outdir_fastq = config['temp_fastq'],
        cores = '4',
        min_read_length=str(options_dict['trim_min_length']),
        min_read_quality=str(options_dict['trim_min_quality']),
        adapter_1 = options_dict['trim_adapter1']

    log:
        log= config['log_dir'] + 'trim_fastq_{sample_id}.log'
    shell:
        '''
        trim_galore --dont_gzip --fastqc --fastqc_args "--outdir {params.outdir_fastqc}" -o {params.outdir_fastq} \
        -a {params.adapter_1}\
        --cores {params.cores} --length {params.min_read_length} --quality {params.min_read_quality} \
        {input.fq1} \
        2> {log.log}
        '''