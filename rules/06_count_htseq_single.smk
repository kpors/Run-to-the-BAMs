rule count_htseq_single:
    # Description
    input:
        bam=config['output_dir'] + date_dir + 'bam_files/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam',
        gtf = config['temp_index'] + gtf_file.replace('.gtf.gz', '.gtf')
    output:
        config['temp_count'] + date_dir + '{sample_id}.tsv'
    conda:
        '../' + config['envs_dir'] + 'htseq.yml'
    #threads:
    resources:
        mem = config['slurm']['htseq']['mem'],
        runtime = config['slurm']['htseq']['runtime']
    params:
        strand_mode=strand_mode,
        align_mode=align_mode
    log:
        log = config['log_dir'] + '{sample_id}_htseq.log'
    shell:
        '''htseq-count --stranded {params.strand_mode} --format bam --mode {params.align_mode} \
        --with-header --order pos --nonunique all\
        --counts_output {output} {input.bam} {input.gtf} \
        2> {log.log}
        '''
