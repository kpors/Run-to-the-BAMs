rule indexing:
    # Description
    input:
        fasta_file = config['temp_index'] + fa_file.replace('.fa.gz', '.fa'),
        gtf_file = config['temp_index'] + gtf_file.replace('.gtf.gz', '.gtf')
    output:
        config['temp_index'] + 'Genome'
    conda:
        '../' + config['envs_dir'] + 'star.yml'
    threads: config['slurm']['star_index']['threads']
    resources:
        mem = config['slurm']['star_index']['mem'],
        runtime = config['slurm']['star_index']['runtime']
    params:
        output_folder = config['temp_index'],
        cores = str(int(config['slurm']['star_index']['threads']) - 1),
        read_length = str(int(options_dict['read_length']) - 1)
    log:
        log = config['log_dir'] + 'indexing.log'
    shell:
        '''
        STAR --runMode genomeGenerate --runThreadN {params.cores} --genomeDir {params.output_folder} \
        --genomeFastaFiles {input.fasta_file} --sjdbGTFfile {input.gtf_file} \
        --sjdbOverhang {params.read_length} \
        2>{log.log}
        '''