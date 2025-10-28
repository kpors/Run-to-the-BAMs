rule rename_fq:
    input:
        fq1 = lambda wc: config['fastq_dir'] + '/' + sample_data_dict[wc.sample_id][0],
        fq2 = lambda wc: config['fastq_dir'] + '/' + sample_data_dict[wc.sample_id][1]
    output:
        fq1 = config['fastq_dir'] + '{sample_id}_1.fastq.gz',
        fq2 = config['fastq_dir'] + '{sample_id}_2.fastq.gz'
    shell:
        """
        ln -s {input.fq1} {output.fq1}
        ln -s {input.fq2} {output.fq2}
        """