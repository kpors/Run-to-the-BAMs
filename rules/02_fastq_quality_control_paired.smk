rule fastq_quality_control_paired:
    # Description
    input:
        input_files = lambda wc:[
            config['fastq_dir'] + sample_data_dict[wc.sample_id][0].replace('fq.gz', 'fastq.gz'),
            config['fastq_dir'] + sample_data_dict[wc.sample_id][1].replace('fq.gz', 'fastq.gz')
        ]
    output:
        config['qc_temp'] + date_dir + '{sample_id}_1_fastqc.zip',
        config['qc_temp'] + date_dir + '{sample_id}_2_fastqc.zip'
    conda:
        '../' + config['envs_dir'] + 'fastqc.yml'
    #threads:
    #resources:
    params:
        out_dir = config['qc_temp'] + date_dir
    log:
        log = config['log_dir'] + '{sample_id}_fastqc.log'
    shell:
        '''fastqc --outdir ./{params.out_dir} ./{input.input_files} \
        2> {log.log}
        '''