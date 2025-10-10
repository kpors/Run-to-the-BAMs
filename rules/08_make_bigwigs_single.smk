rule bigwig_normalised:
    input:
        bam = config['temp_bam'] + date_dir + '{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam',
        bam_index = config['temp_bam'] + date_dir + '{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam.bai',
        size_factors = config['output_dir'] + date_dir + 'counts/deseq2_size_factors.csv'
    output:
        bw = config['output_dir'] + date_dir + 'bigwig/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP_sfNorm_bin5.bw'
    conda:
        '../' + config['envs_dir'] + 'deeptools.yml'
    threads: config['slurm']['bamcoverage']['threads']
    resources:
        mem=config['slurm']['bamcoverage']['mem'],
        runtime=config['slurm']['bamcoverage']['runtime']
    params:
        threads = str(int(config['slurm']['bamcoverage']['threads']) - 1),
        bin_size = 5,
        out_dir= config['output_dir'] + date_dir + 'bigwig/',
        sample_id = '{sample_id}'
    log:
        log=config['log_dir'] + 'bigwig_sfnorm_plus_{sample_id}.log'
    shell:
        '''
        mkdir -p {params.out_dir}
        size_factor=$(awk -F',' -v id="{params.sample_id}" 'NR>1 {{gsub(/"/, "", $1); if ($1 == id) print 1/$2}}' {input.size_factors})
        bamCoverage -p {params.threads} --binSize {params.bin_size} \
            --scaleFactor "$size_factor" \
            -b {input.bam} -o {output.bw} \
            > {log.log}
        '''

rule bigwig_unnormalised:
    input:
        bam = config['temp_bam'] + date_dir + '{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam',
        bam_index= config['temp_bam'] + date_dir + '{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam.bai'
    output:
        bw = config['output_dir'] + date_dir + 'bigwig/{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP_noNorm_bin5.bw'
    conda:
        '../' + config['envs_dir'] + 'deeptools.yml'
    threads: config['slurm']['bamcoverage']['threads']
    resources:
        mem=config['slurm']['bamcoverage']['mem'],
        runtime=config['slurm']['bamcoverage']['runtime']
    params:
        threads = str(int(config['slurm']['bamcoverage']['threads']) - 1),
        bin_size = 5,
        out_dir= config['output_dir'] + date_dir + 'bigwig/',
        sample_id = '{sample_id}'
    log:
        log=config['log_dir'] + 'bigwig_nonorm_minus_{sample_id}.log'
    shell:
        '''
        mkdir -p {params.out_dir}
        bamCoverage -p {params.threads} --binSize {params.bin_size} \
            -b {input.bam} -o {output.bw} \
            > {log.log}
        '''