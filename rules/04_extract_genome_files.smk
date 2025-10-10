

rule download_fasta:
    # Description
    #input:
    output:
        fa = config['input_genome'] + fa_file
    #conda:
    #threads:
    #resources:
    params:
        fa_link = fa_link
    #log:
    shell:
        '''
        wget -O {output.fa} {params.fa_link}
        '''

rule download_gtf:
    # Description
    #input:
    output:
        gtf = config['input_genome'] + gtf_file
    #conda:
    #threads:
    #resources:
    params:
        gtf_link = gtf_link
    #log:
    shell:
        '''
        wget -O {output.gtf} {params.gtf_link}
        '''

rule unzip_fasta:
    # Description
    input:
        fa = config['input_genome'] + fa_file
    output:
        fa = config['temp_index'] + fa_file.replace('.fa.gz', '.fa')
    #conda:
    #threads:
    #resources:
    #params:
    #log:
    shell:
        '''
        gzip -c --decompress --keep {input.fa} > {output.fa}
        '''

rule unzip_gtf:
    # Description
    input:
        gtf = config['input_genome'] + gtf_file
    output:
        gtf = config['temp_index'] + gtf_file.replace('.gtf.gz', '.gtf')
    #conda:
    #threads:
    #resources:
    #params:
    #log:
    shell:
        '''
        gzip -c --decompress --keep {input.gtf} > {output.gtf}
        '''