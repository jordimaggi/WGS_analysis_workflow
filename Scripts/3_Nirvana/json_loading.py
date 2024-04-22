from datetime import datetime
import gzip
import json
import pandas as pd

#Defining JSON data loading function


def json_loading(f,g):

    begin_time = datetime.now()

    genes_df = pd.read_csv(g, sep='\t')
    genes_list = genes_df['gene symbol'].tolist()

    previous_chrom = ''
    fyle = f.split('/')[-1]
    s = fyle.split('_')[0]
    header = ''
    rare_positions = []
    genes = []
    is_header_line = True
    is_position_line = False
    is_gene_line = False
    gene_section_line = '],"genes":['
    end_line = ']}'
    first_time = True
    variants_field = 'variants'
    gnomad_field = 'gnomad'
    freq_threshold = 0.01


    with gzip.open(f, 'rt') as f:
        rare_position_count = 0
        gene_count = 0
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
                    position_dict = json.loads(trim_line.rstrip(','))
                    chrom = position_dict['chromosome']
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
                                for item in transcripts[0]:
                                    if item in genes_list:

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
                                                rare_positions.append(trim_line.rstrip(','))
                                                rare_position_count += 1
                                        elif gnomad_field not in variant_dict:
                                            rare_positions.append(trim_line.rstrip(','))
                                            rare_position_count += 1
                    previous_chrom = position_dict['chromosome']
                if is_gene_line:
                    ## remove the trailing ',' if there is
                    genes.append(trim_line.rstrip(','))
                    gene_count += 1

    chr_ind = {'chr1': [], 'chr2': [], 'chr3': [], 'chr4': [], 'chr5': [], 'chr6': [], 'chr7': [], 'chr8': [], 'chr9': [], 'chr10': [], 'chr11': [], 'chr12': [], \
          'chr13': [], 'chr14': [], 'chr15': [], 'chr16': [], 'chr17': [], 'chr18': [], 'chr19': [], 'chr20': [], 'chr21': [], 'chr22': [], 'chrX': [], 'chrY': [], 'chrM': [], 'chrMT': []}

    for i in range(rare_position_count):
        position_dict = json.loads(rare_positions[i])
        chrom = position_dict['chromosome']
        if chrom not in chr_ind:
            continue
        if chrom != previous_chrom:
            chr_ind[chrom].append(i)
            if i > 0:
                chr_ind[previous_chrom].append(i-1)
        previous_chrom = position_dict['chromosome']
    chr_ind[previous_chrom].append(i)

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

    return chr_ind, rare_positions, genes, samples
