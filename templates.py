import os

###############
# Sample sheet
###############

def extract_sample_info(sample_sheet):
    """Function to extract all information for sample_sheet in input folder,
    and to check if all checkboxes are correctly filled"""
   # Read sample table
    with open(sample_sheet, 'r') as sample_file:
        lines = sample_file.readlines()

        # Collect option inputs and control if check box is used
        # Possible values:
        # 0: No string in checkbox
        # 1 or above: String in checkbox
        dict_keys = ['organism', 'fa_link', 'gtf_link', 'seq_type', 'read_length', 'input_type', 'trim_detection',
                     'trim_adapter1', 'trim_adapter2', 'trim_min_length', 'trim_min_quality',
                     'map_multimapper', 'map_mismatches', 'count_strandedness', 'count_align_mode',
                     'normalisation']
        dict_values = []
        # Numeric options
        line_idx = lines[3].split(';')
        dict_value = line_idx[2]
        input_type = dict_value
        if 1 <= int(dict_value) <= 3:
            dict_values.append(str(dict_value))
        else:
            dict_values.append('0')
        # Organism links
        for line in lines[4:6]:
            line_idx = line.split(';')
            dict_value = line_idx[2]
            if dict_value == '':
                dict_values.append('0')
            else:
                dict_values.append(str(dict_value))
        # Numeric options
        for line in lines[6:10]:
            line_idx = line.split(';')
            dict_value = line_idx[2]
            if 1 <= int(dict_value) <= 500:
                dict_values.append(str(dict_value))
            else:
                dict_values.append('0')
        # Adapter options
        for line in lines[10:12]:
            line_idx = line.split(';')
            dict_value = line_idx[2]
            if dict_value == '':
                dict_values.append('0')
            else:
                dict_values.append(str(dict_value))
        # Numeric options
        for line in lines[12:19]:
            line_idx = line.split(';')
            dict_value = line_idx[2]
            if 1 <= int(dict_value) <= 1000:
                dict_values.append(str(dict_value))
            else:
                dict_values.append('0')

        options_dict = dict(zip(dict_keys, dict_values))

        # Collect data inputs
        sample_data_dict = {}  # For data and input control management
        for line in lines[21:]:
            if options_dict['input_type'] == '1':
                line_idx = line.split(';')[1:]
                sample_id = line_idx[0]
                if sample_id == '':
                    break
                else:
                    sample_data_dict[sample_id] = {}
                # Collect data info for sample_id
                sample_data_dict[sample_id] = [str(line_idx[1]), str(line_idx[2]), str(line_idx[3]),
                                               str(line_idx[4])]
            else:
                line_idx = line.split(';')[1:]
                sample_id = line_idx[0]
                if sample_id == '':
                    break
                else:
                    sample_data_dict[sample_id] = {}
                # Collect data info for sample_id
                if options_dict['seq_type'] == '1':
                    sample_data_dict[sample_id] = [sample_id + '_1.fastq.gz', sample_id + '_2.fastq.gz',
                                                   str(line_idx[3]), str(line_idx[4])]
                else:
                    sample_data_dict[sample_id] = [sample_id + '.fastq.gz', str(line_idx[2]),
                                                   str(line_idx[3]), str(line_idx[4])]
    return options_dict, sample_data_dict




def messages(sample_sheet, fastq_dir):
    """ Function to check for errors when filling the sample info sheet and shout out error messages
    if errors are found - or make a information shout out if everything is ok"""
    options_dict, sample_data_dict = extract_sample_info(sample_sheet)

    # Make dict to collect ready-to-rumble information
    run_check = {}

    # Check for errors in options_dict
    options_error = []
    for key in options_dict.keys():
        if options_dict[key] == 0:
            options_error.append(False)
        else:
            options_error.append(True)

    # Collect options checks and make shout out line
    error_dict_options = {}
    good_dict_options = {}
    if all(options_error):
        run_check['options'] = True
        options_error_messages = 'Na'
        if options_dict['organism'] == '1':
            good_dict_options['organism'] = 'Homo sapiens'
        elif options_dict['organism'] == '2':
            good_dict_options['organism'] = 'Mus musculus'
        elif options_dict['organism'] == '3':
            good_dict_options['organism'] = 'Other species'

        if options_dict['seq_type'] == '1':
            good_dict_options['seq_type'] = 'Paired end'
        elif options_dict['seq_type'] == '2':
            good_dict_options['seq_type'] = 'Single end'

        if options_dict['input_type'] == '1':
            good_dict_options['input_type'] = 'Input folder'
        elif options_dict['input_type'] == '2':
            good_dict_options['input_type'] = 'SRA Archive'

        if options_dict['trim_detection'] == '1':
            good_dict_options['trim_detection'] = 'Auto-detection'
        elif options_dict['trim_detection'] == '2':
            good_dict_options['trim_detection'] = 'Manual'

        if options_dict['count_strandedness'] == '1':
            good_dict_options['count_strandedness'] = 'Reverse'
        elif options_dict['count_strandedness'] == '2':
            good_dict_options['count_strandedness'] = 'Forward'
        elif options_dict['count_strandedness'] == '3':
            good_dict_options['count_strandedness'] = 'No'

        if options_dict['count_align_mode'] == '1':
            good_dict_options['count_align_mode'] = 'Intersection strict'
        elif options_dict['count_align_mode'] == '2':
            good_dict_options['count_align_mode'] = 'Union'
        elif options_dict['count_align_mode'] == '3':
            good_dict_options['count_align_mode'] = 'Intersection-nonempty'

        if options_dict['normalisation'] == '1':
            good_dict_options['normalisation'] = 'TPM'
        elif options_dict['normalisation'] == '2':
            good_dict_options['normalisation'] = 'DeSeq2'
        elif options_dict['normalisation'] == '3':
            good_dict_options['normalisation'] = 'TPM & DeSeq2'
        options_good_messages = '\n' \
        '\n' \
        '###########################################################' + '\n' \
        '# Options for RNA sequencing data processing              #' + '\n' \
        '###########################################################' + '\n' \
        '# Organism:                           ' + good_dict_options['organism'] + '\n' \
        '# Sequencing type:                    ' + good_dict_options['seq_type'] + '\n' \
        '# Input type:                         ' + 'From ' + good_dict_options['input_type'] + '\n' \
        '# Read length:                        ' + str(options_dict['read_length']) + '\n' \
        '\n' \
        '# Trimming parameters:' + '\n' \
        '#   Adapter detection:                ' + good_dict_options['trim_detection'] + '\n' \
        '#   Minimum read length:              ' + options_dict['trim_min_length'] + '\n' \
        '#   Minimum read quality:             ' + options_dict['trim_min_quality'] + '\n' \
        '\n' \
        '# Mapping parameters:' + '\n' \
        '#   Number of multimappers allowed:   ' + options_dict['map_multimapper'] + '\n' \
        '#   Number of mismatches allowed:     ' + options_dict['map_mismatches'] + '\n' \
        '\n' \
        '# Counting and normalisation parameters:' + '\n' \
        '#   Strandedness:                     ' + good_dict_options['count_strandedness'] + '\n' \
        '#   Alignment mode:                   ' + good_dict_options['count_align_mode'] + '\n' \
        '#   Normalisation mode:               ' + good_dict_options['normalisation'] + '\n' \
        '###########################################################' + '\n'


    else:
        run_check['options'] = False
        options_good_messages = 'Na'
        if options_dict['organism'] == '0':
            error_dict_options['organism'] = '- Missing organism information'
        if options_dict['organism'] == '3':
            if options_dict['fa_link'] == '0':
                error_dict_options['fa_link'] = '- Missing FA file link information'
            elif options_dict['gtf_link'] == '0':
                error_dict_options['gtf_link'] = '- Missing GTF file link information'
        if options_dict['seq_type'] == '0':
            error_dict_options['seq_type'] = '- Missing sequencing type information'
        if options_dict['read_length'] == '0':
            error_dict_options['read_length'] = '- Missing read length information'
        if options_dict['input_type'] == '0':
            error_dict_options['input_type'] = '- Missing sequencing data input information'
        if options_dict['trim_detection'] == '0':
            error_dict_options['trim_detection'] = '- Missing trimming adapter detection mode information'
        elif options_dict['trim_detection'] == '2':
            if options_dict['trim_adapter1'] == '0':
                error_dict_options['trim_adapter1'] = '- Missing adapter sequence information'
            if options_dict['seq_type'] == '2':
                if options_dict['trim_adapter2'] == '0':
                    error_dict_options['trim_adapter2'] = '- Missing adapter sequence information'
        if options_dict['trim_min_length'] == '0':
            error_dict_options['trim_min_length'] = '- Missing trimming minimum read length information'
        if options_dict['trim_min_quality'] == '0':
            error_dict_options['trim_min_quality'] = '- Missing trimming minimum read quality information'
        if options_dict['map_multimapper'] == '0':
            error_dict_options['map_multimapper'] = '- Missing multi mapping information'
        if options_dict['map_mismatches'] == '0':
            error_dict_options['map_mismatches'] = '- Missing mapping mismatches information'
        if options_dict['count_strandedness'] == '0':
            error_dict_options['count_strandedness'] = '- Missing counting strandedness information'
        if options_dict['count_align_mode'] == '0':
            error_dict_options['count_align_mode'] = '- Missing counting alignment mode information'
        if options_dict['normalisation'] == '0':
            error_dict_options['normalisation'] = '- Missing counting alignment mode information'
        # Join all errors
        options_error_messages = \
        '###########################################################' + '\n' \
        '# Errors noted - has to be fixed before running pipeline  #' + '\n' \
        '###########################################################' + '\n' \
        '# ERROR in sample sheet information:\n' + ('\n').join(error_dict_options.values()) + '\n' \
        '###########################################################'



    ### Data information ###
    # Check for errors
    # Run only if no errors in options
    if run_check['options']:
        error_samples = []
        good_samples = []
        # Collect input filenames
        input_file_names = os.listdir(fastq_dir)

        # Check fastq files are in input folder if input is needed
        if options_dict['input_type'] == 1:
            for key in sample_data_dict.keys():
                if str(sample_data_dict[key][1]):
                    if sample_data_dict[key][1] in input_file_names:
                        if options_dict['seq_type'] == '1':
                            if str(sample_data_dict[key][2]):
                                if sample_data_dict[key][2] in input_file_names:
                                    good_samples.append(key)
                                else:
                                    error_samples.append(key)
                            else:
                                error_samples.append(key)
                        else:
                            good_samples.append(key)
                    else:
                        error_samples.append(key)
                else:
                    error_samples.append(key)

        # Check design factors are noted if needed
        if options_dict['normalisation'] in range(2,3):
            for key in sample_data_dict.keys():
                if str(sample_data_dict[key][3]):
                    if str(sample_data_dict[key][4]):
                        good_samples.append(key)
                    else:
                        error_samples.append(key)
                else:
                    error_samples.append(key)

        # Collect samples checks and make shout out line
        if len(error_samples) == 0:
            if len(good_samples) == len(sample_data_dict.keys()):
                run_check['samples'] = True
            else:
                run_check['samples'] = False
        else:
            run_check['samples'] = False

        if run_check['samples']:
            sample_error_messages = 'Na'
        else:
            sample_error_messages = \
            '###########################################################' + '\n' \
            '# Errors noted - has to be fixed before running pipeline  #' + '\n' \
            '###########################################################' + '\n' \
            '# Missing FASTQ file or DeSeq2 design factor information ' + '\n' \
            '# in for following samples:' + '\n' \
            + ('\n').join(list(set(error_samples))) + '\n' \
            '###########################################################'


    # Update run check bool
    run_check_bool = all(run_check)

    return run_check_bool, sample_data_dict, options_dict, sample_error_messages, options_error_messages, options_good_messages


fastq_dir ='./input/fastq/'
sample_sheet = './input/Sample_sheet.csv'

run_check_bool, sample_data_dict, options_dict, sample_error_messages, options_error_messages, options_good_messages = messages(sample_sheet, fastq_dir)

# print(sample_ids)
# print(sra_dict)