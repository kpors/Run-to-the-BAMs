rule merge_counts:
    # Description
    input:
        count_files = expand(config['temp_count'] + date_dir + '{sample_id}.tsv', sample_id=sample_ids),
        gtf_file = config['temp_index'] + gtf_file.replace('.gtf.gz','.gtf')
    output:
        config['output_dir'] + date_dir + 'counts/counts_raw.csv'
    conda:
        '../' + config['envs_dir'] + 'r_normalise.yml'
    #threads:
    resources:
        mem=config['slurm']['normalise']['mem'],
        runtime=config['slurm']['normalise']['runtime']
    params:
        sample_names = (',').join(sample_ids),
        out_dir = config['output_dir'] + date_dir + 'counts/',
        countfile_dir = config['temp_count'] + date_dir
    log:
        log = config['log_dir'] + 'merge_counts.log'
    shell:
        '''
        mkdir -p {params.out_dir}
        Rscript ./scripts/merge_counts.R {params.countfile_dir} {params.sample_names} {input.gtf_file} {output} \
        2> {log.log}
        '''