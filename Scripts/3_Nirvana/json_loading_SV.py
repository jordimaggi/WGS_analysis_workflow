from datetime import datetime
import gzip
import json
import pandas as pd

#Defining JSON data loading function


def json_loading(f,g,s,fold):

    begin_time = datetime.now()

    genes_df = pd.read_csv(g, sep='\t')
    genes_list = genes_df['gene symbol'].tolist()

    #Creating a list containing the IDs from the VCF files
    file = '{}/{}/SV_calling/{}_SV_merged.norm.vcf.gz'.format(fold, s, s)
    with gzip.open(file) as file_in:
        vcf = pd.read_csv(file_in, sep='\t', header=326)
    IDs_raw = vcf['ID'].tolist()
    del vcf

    previous_chrom = ''
    fyle = f.split('/')[-1]
    s = fyle.split('_')[0]
    header = ''
    rare_positions = []
    genes = []
    IDs = []
    is_header_line = True
    is_position_line = False
    is_gene_line = False
    gene_section_line = '],"genes":['
    end_line = ']}'
    first_time = True
    variants_field = 'variants'
    gnomad_field = 'gnomad'
    freq_threshold = 0.01
    tot_variants_count = 0


    with gzip.open(f, 'rt') as f:
        rare_position_count = 0
        gene_count = 0
        count = 0
        for line in f:

            trim_line = line.strip()
            if is_header_line:
                ## only keep the "header" field content from the line
                header = trim_line[10:-14]
                is_header_line = False
                is_position_line = True
                continue
            if trim_line == gene_section_line:
                is_gene_line = True
                is_position_line = False
                continue
            elif trim_line == end_line:
                break
            else:
                if is_position_line:
                    try:
                        position_dict = json.loads(trim_line.rstrip(','))
                        chrom = position_dict['chromosome']
                        if chrom == 'chr2' and position_dict['position'] > 33100000 and position_dict['position'] < 33150000 and '33141' in str(position_dict['altAlleles']):
                            continue
                        if chrom != previous_chrom:
                            print('Collecting rare variants on {}'.format(chrom))
                        if variants_field in position_dict:
                            for variant_dict in position_dict[variants_field]:
                                if 'transcripts' in variant_dict:
                                    transcripts = [[]]
                                    for trans in variant_dict['transcripts']:
                                        if trans['source'] == 'Ensembl':
                                            continue
                                            #Keeping only RefSeq entries
                                        else:
                                            transcripts[0].append(trans['hgnc'])
                                    transcripts[0] = list(dict.fromkeys(transcripts[0]))
                                    added = False
                                    for item in transcripts[0]:
                                        if item in genes_list and added == False:
                                            if gnomad_field in variant_dict:
                                                if 'allAf' not in variant_dict[gnomad_field]:
                                                    freq = 0
                                                else:
                                                    freq = variant_dict[gnomad_field]['allAf']
                                                if freq == '':
                                                    freq = 0
                                                if freq > freq_threshold:
                                                    continue

                                                else:
                                                    if trim_line.rstrip(',') not in rare_positions:
                                                        rare_positions.append(trim_line.rstrip(','))
                                                        rare_position_count += 1
                                                        IDs.append(IDs_raw[tot_variants_count])
                                                        added = True
                                            elif gnomad_field not in variant_dict:
                                                if trim_line.rstrip(',') not in rare_positions:
                                                    rare_positions.append(trim_line.rstrip(','))
                                                    IDs.append(IDs_raw[tot_variants_count])
                                                    rare_position_count += 1
                                                    added = True
                        previous_chrom = position_dict['chromosome']
                        tot_variants_count += 1
                    except:
                        count += 1
                if is_gene_line:
                    ## remove the trailing ',' if there is
                    genes.append(trim_line.rstrip(','))
                    gene_count += 1


    print(' ')
    print('Lines with errors in the annotations that could not be processed: {}'.format(count))
    print(' ')
    print ('Annotated using the following sources:')
    lyst = header.split('[{"name":')[1].split('},{"name":')
    for el in lyst:
        print(el)

    samples = header.split('[{"name":')[1].split('},{"name":')[-1].split('],"samples":[')[1].split(']}')[0].split(', ')[0].split('","')
    samples[0] = samples[0][1:]
    samples[-1] = samples[-1][:-1]
    print('Annotated files contain variants in samples: ', samples)
    print(' ')
    print('number of rare variants (< 1% gnomAD all):', rare_position_count)
    print('number of genes:', gene_count)
    print(' ')
    print('Parsing through the {} file for sample {} was finished in {}.'.format(fyle, s, datetime.now() - begin_time))
    print(' ')

    return rare_positions, genes, samples, IDs
