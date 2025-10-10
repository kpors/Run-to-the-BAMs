rule fastq_download_single_end:
    # Description
    #input:
    output:
        fq1 = config['fastq_dir'] + '{sample_id}.fastq.gz'
    conda: '../' + config['envs_dir'] + 'sra_tools.yml'
    threads: config['slurm']['sra_tools']['threads']
    resources:
        mem=config['slurm']['sra_tools']['mem'],
        runtime=config['slurm']['sra_tools']['runtime']
    params:
        sra_id = '{sample_id}',
        outdir = config['fastq_dir'],
        fq1= config['fastq_dir'] + '{sample_id}.fastq',
        threads = str(int(config['slurm']['sra_tools']['threads']) - 1)
    log:
        log=config['log_dir'] + '{sample_id}_fasterq_dump.log'
    shell:
        '''
        fasterq-dump {params.sra_id} --outdir {params.outdir} -e {params.threads} \
        > {log.log}
        pigz -p {params.threads} {params.fq1}
        '''