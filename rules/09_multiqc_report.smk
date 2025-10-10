rule multiqc_report:
    # Description
    input:
        input_files=expand(
            config['temp_bam'] + date_dir + '{sample_id}_' + str(genome_id) + '_' + str(multimap_id) + '_noDDUP.bam',
            sample_id=sample_ids)
    output:
        output_files = config['output_dir'] + date_dir + 'multiqc_report.html'
    conda:
        '../' + config['envs_dir'] + 'fastqc.yml'
    #threads:
    #resources:
    params:
        fastqc_dir=config['qc_temp'] + date_dir,
        bam_dir=config['temp_bam'] + date_dir,
        out_dir = config['output_dir'] + date_dir
    log:
        log = config['log_dir'] + 'multiqc_report.log'
    shell:
        '''multiqc --outdir ./{params.out_dir} ./{params.fastqc_dir} ./{params.bam_dir} \
        2> {log.log}
        '''