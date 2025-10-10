rule normalise_deseq2:
    # Description
    input:
        count_file = config['output_dir'] + date_dir + 'counts/counts_raw.csv'
    output:
        counts = config['output_dir'] + date_dir + 'counts/counts_normalised_deseq2.csv',
        sf = config['output_dir'] + date_dir + 'counts/deseq2_size_factors.csv'
    conda:
        '../' + config['envs_dir'] + 'r_normalise.yml'
    #threads:
    #resources:
    params:
        out_dir = config['output_dir'] + date_dir + 'counts/'
    log:
        log = config['log_dir'] + 'normalise_deseq2.log'
    shell:
        '''
        Rscript ./scripts/normalise_deseq2.R {input.count_file} {params.out_dir} \
        2> {log.log}
        '''

rule normalise_tpm:
    # Description
    input:
        count_file = config['output_dir'] + date_dir + 'counts/counts_raw.csv',
        gtf_file = config['temp_index'] + gtf_file.replace('.gtf.gz','.gtf')
    output:
        counts = config['output_dir'] + date_dir + 'counts/counts_normalised_tpm.csv'
    conda:
        '../' + config['envs_dir'] + 'r_normalise.yml'
    #threads:
    resources:
        mem=config['slurm']['normalise']['mem'],
        runtime=config['slurm']['normalise']['runtime']
    params:
        out_dir = config['output_dir'] + date_dir + 'counts/'
    log:
        log = config['log_dir'] + 'normalise_tpm.log'
    shell:
        '''
        Rscript ./scripts/normalise_tpm.R {input.count_file} {input.gtf_file} {output.counts} \
        2> {log.log}
        '''