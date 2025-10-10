
rule align_star_single:
    # Description
    input:
        fq1 = config['temp_fastq'] + '{sample_id}_trimmed.fq',
        index = config['temp_index'] + 'Genome'
    output:
        config['temp_bam'] + date_dir + '{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam'
    conda: '../' + config['envs_dir'] + 'star.yml'
    threads: config['slurm']['star_align']['threads']
    resources:
        mem=config['slurm']['star_align']['mem'],
        runtime=config['slurm']['star_align']['runtime']
    params:
        cores = str(int(config['slurm']['star_align']['threads']) - 1),
        index_base = config['temp_index'],
        mismatches_max=str(options_dict['map_mismatches']),
        n_multimappers=str(options_dict['map_multimapper']),
        gtf_file = config['temp_index'] + gtf_file.replace('.gtf.gz', '.gtf'),
        read_length = str(int(options_dict['read_length']) - 1),
        output_prefix = config['temp_bam'] + date_dir + '{sample_id}_',
        output_name = config['temp_bam'] + date_dir + '{sample_id}_Aligned.sortedByCoord.out.bam'
    log:
        log= config['log_dir'] + 'star_align_{sample_id}.log'
    shell:
        '''
        STAR --runMode alignReads --runThreadN {params.cores} --genomeDir {params.index_base} \
        --outFileNamePrefix {params.output_prefix} --sjdbGTFfile {params.gtf_file} \
        --winAnchorMultimapNmax {params.n_multimappers} --outFilterMultimapNmax {params.n_multimappers} \
        --outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonical \
        --outFilterMismatchNmax {params.mismatches_max} --sjdbOverhang {params.read_length} --outReadsUnmapped Fastx \
        --readFilesIn {input.fq1} --quantMode GeneCounts \
        2> {log.log}
        mv {params.output_name} {output}
        '''


rule index_bam:
    # Description
    input:
        bam = config['temp_bam'] + date_dir + '{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam'
    output:
        index=config['temp_bam'] + date_dir + '{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam.bai'
    conda: '../' + config['envs_dir'] + 'samtools.yml'
    #threads:
    #resources:
    #params:
    #log:
    shell:
        '''
        samtools index {input.bam}
        '''