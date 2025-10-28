rule save_pipeline:
    # Description
    #input:
    output:
        output_files = config['output_dir'] + date_dir + 'pipeline_parameters/templates.py'
    #conda:
    #threads:
    #resources:
    params:
        rules_dir='rules/',
        envs_dir=config['envs_dir'],
        config='config/config.yaml',
        snakefile='Snakefile',
        templatesfile='templates.py',
        samplefile=config['sample_table'],
        scripts=config['scripts_dir'],
        out_dir = config['output_dir'] + date_dir + 'pipeline_parameters/'
    #log:
    shell:
        '''
        mkdir -p {params.out_dir}
        cp {params.templatesfile} {params.out_dir}
        cp -r {params.rules_dir} {params.out_dir}
        cp -r {params.envs_dir} {params.out_dir}
        cp {params.config} {params.out_dir}
        cp {params.snakefile} {params.out_dir}
        cp {params.samplefile} {params.out_dir}
        cp -r {params.scripts} {params.out_dir}
        '''