from json_loading import *
from datetime import datetime
import pandas as pd
import os
import json
import gzip
import pandas as pd
import math
import sys

sample = snakemake.wildcards.sample
folder = snakemake.wildcards.folder
sys.stdout = open("{}/{}/Console_outputs/Nirvana_JsonDecoding.{}.out".format(folder, sample, sample), 'w')

start = datetime.now()
print(start)
print("Starting JSON files decoding and filtering command for sample {}".format(sample))
print(' ')

#Loading JSON data and creating chr indeces dictionaries and variants and genes data with the json_loading function
chr_ind_HC, HC_positions, HC_genes = json_loading(snakemake.input.hc)
chr_ind_DV, DV_positions, DV_genes = json_loading(snakemake.input.dv)
print('Haplotypecaller positions for all chromosomes:', chr_ind_HC)
print('Deepvariant positions for all chromosomes:', chr_ind_DV)

#Source: https://www.genscript.com/tools/codon-frequency-table
#Expression Host Organism: Human (Genetic code: Standard)
codon_usage = {'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13, 'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, \
              'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.28, 'TAG': 0.20, 'TGT': 0.45, 'TGC': 0.55, 'TGA': 0.52,'TGG': 1.00, \
              'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.41, 'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11, \
              'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75, 'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, \
              'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1.00, 'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12, \
              'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58, 'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.20, 'AGG': 0.20, \
              'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47, 'GCT': 0.26, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11, \
              'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25}

#Defining recording dictionary for inheritance modes
inher_conv = {'Autosomal dominant':'AD', 'Autosomal recessive':'AR', 'X-linked':'XL', 'X-linked recessive':'XLR', \
              'X-linked dominant':'XLD', 'Mitochondrial':'Mito', 'Digenic recessive': 'DigenicR', 'Digenic dominant': 'DigenicD', \
              'Somatic mutation': 'Somatic', 'Isolated cases': 'Isolated', 'Multifactorial': 'Multifactorial'}

#Creating a list of gene identifiers from the genes of interest list for eye disorders
genes_file = snakemake.input.genes_list
genes_df = pd.read_csv(genes_file, sep='\t')
genes_list = genes_df['gene symbol'].tolist()

#Creating the omim and gnmoAD_scores dictionaries
omim_genes = {}
gnomAD_scores = {}
#Merge HC and DV genes into a single list
genes = HC_genes + DV_genes
genes = list(dict.fromkeys(genes))
for i4 in range(len(genes)):
    gene_dict = json.loads(genes[i4])
    if 'omim' in gene_dict:
        entries = []
        if 'phenotypes' in gene_dict['omim'][0]:
            for el in gene_dict['omim'][0]['phenotypes']:
                text = ''
                if 'inheritances' in el:
                    for i3 in range(len(el['inheritances'])):
                        if el['inheritances'][i3] in inher_conv:
                            text += inher_conv[el['inheritances'][i3]]
                        else:
                            text += el['inheritances'][i3]
                        if (i3+1) != len(el['inheritances']):
                            text += '/'

                    entries.append('{} - {}'.format(el['phenotype'], text))
                else:
                    entries.append('{} - Unspecified'.format(el['phenotype']))
            omim_genes[gene_dict['name']] = entries
        elif 'phenotypes' not in gene_dict['omim'][0]:
            omim_genes[gene_dict['name']] = ['No associated phenotypes']
    elif 'omim' not in gene_dict:
        omim_genes[gene_dict['name']] = ['No entries']
    if 'gnomAD' in gene_dict:
        gnomAD_scores[gene_dict['name']] = gene_dict['gnomAD']

#Defining the variables that will hold the variants data before writing them to files
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', \
          'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

variants_field = 'variants'
gnomad_field = 'gnomad'
freq_data = {'chromosome': [], 'position': [], 'refAllele': [], 'altAllele': [], 'variant_id': [], 'caller': [], 'quality': [], 'fisherStrandBias': [], 'mappingQuality': [], 'alleleDepths': [], 'variantFrequencies': [], \
             'totalDepth': [], 'genotype': [], 'genotypeQuality': [], 'variantType': [],  'hgvsg': [], \
             'gnomAD_all': [], 'gnomAD_NFE': [], 'gnomAD_max': [], 'gnomAD_max_ethnicity': [],'gnomAD_controls': [], 'gnomAD_homCount': [], \
             'gene(s)': [], 'transcript_consequence': [], 'hgsvc': [], 'hgsvp': [], 'codons': [], 'codonsFraction': [], 'codonsUsageBias': [], 'phyloP': [], \
             'RegBaseV1_RegPhred': [], 'RegBaseV1_PatPhred': [], 'CADD_raw': [], 'CADD_PHRED': [], 'DANN': [], 'revel': [], 'primateAI': [], 'LINSIGHT': [], 'FATHMM-XF': [], 'ncER_perc': [], 'NCBoost_score': [], 'NCBoost_chr_rank_perc': [], \
             'spliceAI_accGainScore': [], 'spliceAI_accGainDist': [], 'spliceAI_accLossScore': [], 'spliceAI_accLossDist': [], 'spliceAI_donGainScore': [], 'spliceAI_donGainDist': [], 'spliceAI_donLossScore': [], 'spliceAI_donLossDist': [], \
             'polyPhen': [], 'sift': [], 'OMIM': [], 'clinvar_[Phen - Signif]': [], 'gnomAD_gene_pLi': [], 'gnomAD_gene_pRec': [], 'gnomAD_gene_pNull': [], 'gnomAD_gene_synZ': [], 'gnomAD_gene_misZ': [], 'gnomAD_gene_loeuf': [], \
             'regulatoryRegionID': [], 'regulatoryRegionType': [], 'inLowComplexityRegion': [], 'dbSNP': [], 'MITOMAP_diseases': [], 'MITOMAP_clinSignificance': []}


begin_time = datetime.now()

#Checking if filtered table already exists; if so, old file needs to be deleted or renamed prior to repeating the procedure
if os.path.exists('{}/{}/{}_Nirvana_combined_rare_variants.tsv'.format(folder, sample, sample)):
    print('Files {}_Nirvana_combined_rare_variants.tsv and {}_Nirvana_combined_rare_variants_genes_list.tsv already exists. Delete or rename file before proceeding.'.format(sample, sample))
    exists = True

#If it does not exists, it creates it with the freq_data dictionary keys as column labels
else:
    #Creating file with header
    freq_df = pd.DataFrame({key:pd.Series(value) for key, value in freq_data.items()})
    freq_df.to_csv('{}/{}_Nirvana_combined_rare_variants.tsv'.format(folder, sample, sample), sep='\t', index=False)
    freq_df.to_csv('{}/{}_Nirvana_combined_rare_variants_genes_list.tsv'.format(folder, sample, sample), sep='\t', index=False)
    exists = False


count = 0
#The script processes one chromosome at a time to save RAM space and simplify merging process
for chrom in chroms:
    if exists == True:
        break

    HC_range = chr_ind_HC[chrom]
    DV_range = chr_ind_DV[chrom]
    positions = {}

    if len(HC_range) == 0:
        print('Ranges for {} not present'.format(chrom))
        continue

    #Create a list of positions in the DV VCF to compare when creating the final dictionary to populate the 'caller' field
    for e in range(HC_range[0], HC_range[1]):
        pos = json.loads(HC_positions[e])['position']
        ref = json.loads(HC_positions[e])['refAllele']
        alt = json.loads(HC_positions[e])['altAlleles']
        positions[pos] = ['HC', [e], ref, alt]
    if len(DV_range) > 1:
        for e2 in range(DV_range[0], DV_range[1]):
            pos = json.loads(DV_positions[e2])['position']
            ref = json.loads(DV_positions[e2])['refAllele']
            alt = json.loads(DV_positions[e2])['altAlleles']
            if pos not in positions:
                positions[pos] = ['DV', [e2], ref, alt]
            else:
                if positions[pos][2] == ref and positions[pos][3] == alt:
                    positions[pos][0] = 'HC, DV'
                    positions[pos][1].append(e2)

    positions_l = list(positions.keys())
    positions_l.sort()

    freq_data = {'chromosome': [], 'position': [], 'refAllele': [], 'altAllele': [], 'variant_id': [], 'caller': [], 'quality': [], 'fisherStrandBias': [], 'mappingQuality': [], 'alleleDepths': [], 'variantFrequencies': [], \
             'totalDepth': [], 'genotype': [], 'genotypeQuality': [], 'variantType': [],  'hgvsg': [], \
             'gnomAD_all': [], 'gnomAD_NFE': [], 'gnomAD_max': [], 'gnomAD_max_ethnicity': [],'gnomAD_controls': [], 'gnomAD_homCount': [], \
             'gene(s)': [], 'transcript_consequence': [], 'hgsvc': [], 'hgsvp': [], 'codons': [], 'codonsFraction': [], 'codonsUsageBias': [], 'phyloP': [], \
             'RegBaseV1_RegPhred': [], 'RegBaseV1_PatPhred': [], 'CADD_raw': [], 'CADD_PHRED': [], 'DANN': [], 'revel': [], 'primateAI': [], 'LINSIGHT': [], 'FATHMM-XF': [], 'ncER_perc': [], 'NCBoost_score': [], 'NCBoost_chr_rank_perc': [], \
             'spliceAI_accGainScore': [], 'spliceAI_accGainDist': [], 'spliceAI_accLossScore': [], 'spliceAI_accLossDist': [], 'spliceAI_donGainScore': [], 'spliceAI_donGainDist': [], 'spliceAI_donLossScore': [], 'spliceAI_donLossDist': [], \
             'polyPhen': [], 'sift': [], 'OMIM': [], 'clinvar_[Phen - Signif]': [], 'gnomAD_gene_pLi': [], 'gnomAD_gene_pRec': [], 'gnomAD_gene_pNull': [], 'gnomAD_gene_synZ': [], 'gnomAD_gene_misZ': [], 'gnomAD_gene_loeuf': [], \
             'regulatoryRegionID': [], 'regulatoryRegionType': [], 'inLowComplexityRegion': [], 'dbSNP': [], 'MITOMAP_diseases': [], 'MITOMAP_clinSignificance': []}

    in_list = []
    for i in range(len(positions_l)):
        count += 1
        if count % 10000 == 0:
            print('Processed {} variants in {}.'.format(count, datetime.now() - begin_time))
        if positions[positions_l[i]][0] != 'DV':
            position_dict = json.loads(HC_positions[positions[positions_l[i]][1][0]])
        elif positions[positions_l[i]][0] == 'DV':
            position_dict = json.loads(DV_positions[positions[positions_l[i]][1][0]])
        if variants_field in position_dict:
            for variant_dict in position_dict[variants_field]:
                if gnomad_field in variant_dict:
                    freq = variant_dict[gnomad_field]['allAf']

                #Cleaning the 'transcripts' column
                if 'transcripts' in variant_dict:
                    transcripts = [[],[],[],[],[], []]
                    codons = []
                    for trans in variant_dict['transcripts']:
                        if trans['source'] == 'Ensembl':
                            continue
                            #Keeping only RefSeq entries
                        else:

                            transcripts[0].append(trans['hgnc'])
                            transcripts[1].append(', '.join([str(item) for item in trans['consequence']]))
                            if 'hgvsc' in trans:
                                transcripts[2].append(trans['hgvsc'])
                            if 'hgvsp' in trans:
                                transcripts[3].append(trans['hgvsp'])
                            if 'polyPhenPrediction' in trans:
                                transcripts[4].append(trans['polyPhenPrediction'])
                            if 'siftPrediction' in trans:
                                transcripts[5].append(trans['siftPrediction'])
                            if 'codons' in trans:
                                codons.append(trans['codons'])
                    transcripts[0] = list(dict.fromkeys(transcripts[0]))
                    codons = list(dict.fromkeys(codons))


                #Cleaning the 'clinvar' column
                clinvar = []
                if 'clinvar' in variant_dict:
                    for entry in variant_dict['clinvar']:
                        if 'phenotypes' not in entry:
                            clinvar.append('[No phenotype provided' + ' - ' + ', '.join([str(item) for item in entry['significance']]) + ']')
                        else:
                            clinvar.append('[' + ', '.join([str(item) for item in entry['significance']]) + ' - ' + ', '.join([str(item) for item in entry['significance']]) + ']')
                clinvar = list(dict.fromkeys(clinvar))

                #Cleaning the spliceAI data
                spliceai = []
                if 'spliceAI'  in variant_dict:
                    if 'acceptorGainScore' not in variant_dict['spliceAI'][0]:
                        spliceai.append('')
                        spliceai.append('')
                    else:
                        spliceai.append(variant_dict['spliceAI'][0]['acceptorGainScore'])
                        spliceai.append(variant_dict['spliceAI'][0]['acceptorGainDistance'])
                    if 'acceptorLossScore' not in variant_dict['spliceAI'][0]:
                        spliceai.append('')
                        spliceai.append('')
                    else:
                        spliceai.append(variant_dict['spliceAI'][0]['acceptorLossScore'])
                        spliceai.append(variant_dict['spliceAI'][0]['acceptorLossDistance'])
                    if 'donorGainScore' not in variant_dict['spliceAI'][0]:
                        spliceai.append('')
                        spliceai.append('')
                    else:
                        spliceai.append(variant_dict['spliceAI'][0]['donorGainScore'])
                        spliceai.append(variant_dict['spliceAI'][0]['donorGainDistance'])
                    if 'donorLossScore' not in variant_dict['spliceAI'][0]:
                        spliceai.append('')
                        spliceai.append('')
                    else:
                        spliceai.append(variant_dict['spliceAI'][0]['donorLossScore'])
                        spliceai.append(variant_dict['spliceAI'][0]['donorLossDistance'])

                sample_dict = position_dict['samples'][0]
                gnomad_ethn = {'afrAf':'AFR', 'amrAf':'AMR', 'easAf':'EAS', 'finAf':'FIN', 'nfeAf':'NFE', 'asjAf':'ASJ', 'othAf':'OTH'}
                gnomad_max = [0, '']
                if gnomad_field in variant_dict:
                    for key in variant_dict[gnomad_field]:
                        if key in gnomad_ethn:
                            if key=='nfeAf':
                                gnomad_NFE = variant_dict[gnomad_field][key]
                            if variant_dict[gnomad_field][key] > gnomad_max[0]:
                                gnomad_max[0] = variant_dict[gnomad_field][key]
                                gnomad_max[1] = gnomad_ethn[key]
                        else:
                            continue

                freq_data['chromosome'].append(position_dict['chromosome'])
                freq_data['position'].append(position_dict['position'])
                freq_data['refAllele'].append(position_dict['refAllele'])
                freq_data['altAllele'].append(', '.join([str(item) for item in position_dict['altAlleles']]))
                freq_data['variant_id'].append(variant_dict['vid'])
                freq_data['caller'].append(positions[positions_l[i]][0])
                if len(positions[positions_l[i]][0]) == 2:
                    freq_data['quality'].append(position_dict['quality'])
                    if 'alleleDepths' in sample_dict:
                        freq_data['alleleDepths'].append(sample_dict['alleleDepths'])
                        freq_data['variantFrequencies'].append(', '.join([str(item) for item in sample_dict['variantFrequencies']]))
                        freq_data['totalDepth'].append(sample_dict['totalDepth'])
                    else:
                        freq_data['alleleDepths'].append('')
                    freq_data['genotype'].append(sample_dict['genotype'])
                    freq_data['genotypeQuality'].append(sample_dict['genotypeQuality'])
                else:
                    DV_qual = json.loads(DV_positions[positions[positions_l[i]][1][1]])['quality']
                    freq_data['quality'].append(str(position_dict['quality']) + ', ' + str(DV_qual))
                    if 'alleleDepths' in json.loads(DV_positions[positions[positions_l[i]][1][1]])['samples'][0]:
                        DV_allDepths = json.loads(DV_positions[positions[positions_l[i]][1][1]])['samples'][0]['alleleDepths']
                        freq_data['alleleDepths'].append('[' + ', '.join([str(item) for item in sample_dict['alleleDepths']]) + '], [' + ', '.join([str(item) for item in DV_allDepths]) + ']')
                        DV_variantFreq = json.loads(DV_positions[positions[positions_l[i]][1][1]])['samples'][0]['variantFrequencies']
                        freq_data['variantFrequencies'].append(', '.join([str(item) for item in sample_dict['variantFrequencies']]) + ', ' + ', '.join([str(item) for item in DV_variantFreq]))
                        DV_totalDep = json.loads(DV_positions[positions[positions_l[i]][1][1]])['samples'][0]['totalDepth']
                        freq_data['totalDepth'].append(str(sample_dict['totalDepth']) + ', ' + str(DV_totalDep))
                    DV_geno = json.loads(DV_positions[positions[positions_l[i]][1][1]])['samples'][0]['genotype']
                    freq_data['genotype'].append(sample_dict['genotype'] + ', ' + DV_geno)
                    DV_genoqual = json.loads(DV_positions[positions[positions_l[i]][1][1]])['samples'][0]['genotypeQuality']
                    freq_data['genotypeQuality'].append(str(sample_dict['genotypeQuality']) + ', ' + str(DV_genoqual))
                if 'fisherStrandBias' in position_dict:
                    freq_data['fisherStrandBias'].append(position_dict['fisherStrandBias'])
                else:
                    freq_data['fisherStrandBias'].append('')
                if 'mappingQuality' in position_dict:
                    freq_data['mappingQuality'].append(position_dict['mappingQuality'])
                else:
                    freq_data['mappingQuality'].append('')
                freq_data['variantType'].append(variant_dict['variantType'])
                if 'hgvsg' in variant_dict:
                    freq_data['hgvsg'].append(variant_dict['hgvsg'])
                else:
                    freq_data['hgvsg'].append('')
                if gnomad_field in variant_dict:
                    freq_data['gnomAD_all'].append(freq)
                    freq_data['gnomAD_max'].append(gnomad_max[0])
                    freq_data['gnomAD_max_ethnicity'].append(gnomad_max[1])
                    try:
                        freq_data['gnomAD_NFE'].append(gnomad_NFE)
                    except:
                        freq_data['gnomAD_NFE'].append('')
                    freq_data['gnomAD_homCount'].append(variant_dict[gnomad_field]['allHc'])
                    if 'controlsAllAf' in variant_dict[gnomad_field]:
                        freq_data['gnomAD_controls'].append(variant_dict[gnomad_field]['controlsAllAf'])
                    else:
                        freq_data['gnomAD_controls'].append('')
                else:
                    freq_data['gnomAD_all'].append('')
                    freq_data['gnomAD_max'].append(gnomad_max[0])
                    freq_data['gnomAD_max_ethnicity'].append('')
                    freq_data['gnomAD_NFE'].append('')
                    freq_data['gnomAD_homCount'].append('')
                    freq_data['gnomAD_controls'].append('')
                if 'transcripts' in variant_dict:
                    freq_data['gene(s)'].append(', '.join([str(item) for item in transcripts[0]]))
                    for item in transcripts[0]:
                        OMIM = ''
                        pLi = []
                        pRec = []
                        pNull = []
                        synZ = []
                        misZ = []
                        loeuf = []
                        if item in genes_list:
                            in_list.append(i)
                        if item in omim_genes:
                            if len(OMIM) > 1:
                                OMIM += '/'
                            OMIM += ', '.join([str(item2) for item2 in omim_genes[item]])
                        else:
                            if len(OMIM) > 1:
                                continue
                            else:
                                OMIM += 'No entries'
                        if item in gnomAD_scores:
                            try:
                                pLi.append(gnomAD_scores[item]['pLi'])
                                pRec.append(gnomAD_scores[item]['pRec'])
                                pNull.append(gnomAD_scores[item]['pNull'])
                                loeuf.append(gnomAD_scores[item]['loeuf'])
                            except:
                                pLi.append('NA')
                                pRec.append('NA')
                                pNull.append('NA')
                                loeuf.append('NA')
                            try:
                                synZ.append(gnomAD_scores[item]['synZ'])
                                misZ.append(gnomAD_scores[item]['misZ'])
                            except:
                                synZ.append('NA')
                                misZ.append('NA')
                                print(gnomAD_scores[item])
                    freq_data['OMIM'].append(OMIM)
                    freq_data['gnomAD_gene_pLi'].append('/'.join([str(item) for item in pLi]))
                    freq_data['gnomAD_gene_pRec'].append('/'.join([str(item) for item in pRec]))
                    freq_data['gnomAD_gene_pNull'].append('/'.join([str(item) for item in pNull]))
                    freq_data['gnomAD_gene_synZ'].append('/'.join([str(item) for item in synZ]))
                    freq_data['gnomAD_gene_misZ'].append('/'.join([str(item) for item in misZ]))
                    freq_data['gnomAD_gene_loeuf'].append('/'.join([str(item) for item in loeuf]))
                    freq_data['transcript_consequence'].append(', '.join([str(item) for item in transcripts[1]]))
                    freq_data['hgsvc'].append(', '.join([str(item) for item in transcripts[2]]))
                    freq_data['hgsvp'].append(', '.join([str(item) for item in transcripts[3]]))
                    freq_data['polyPhen'].append(', '.join([str(item) for item in transcripts[4]]))
                    freq_data['sift'].append(', '.join([str(item) for item in transcripts[5]]))
                    if len(codons) > 0:
                        freq_data['codons'].append(', '.join([str(item) for item in codons]))
                        usages = []
                        biases = []
                        for el in codons:
                            ref_codon = el.split('/')[0].upper()
                            alt_codon = el.split('/')[1].upper()
                            try:
                                ref_usage = codon_usage[ref_codon]
                                alt_usage = codon_usage[alt_codon]
                                usages.append(str(ref_usage) + '/' + str(alt_usage))
                                biases.append("{:.2f}".format(ref_usage - alt_usage))
                            except:
                                usages.append('NA')
                                biases.append('NA')
                        freq_data['codonsFraction'].append(', '.join([str(item) for item in usages]))
                        freq_data['codonsUsageBias'].append(', '.join([str(item) for item in biases]))
                    else:
                        freq_data['codons'].append('')
                        freq_data['codonsFraction'].append('')
                        freq_data['codonsUsageBias'].append('')
                else:
                    freq_data['gene(s)'].append('')
                    freq_data['transcript_consequence'].append('')
                    freq_data['hgsvc'].append('')
                    freq_data['hgsvp'].append('')
                    freq_data['polyPhen'].append('')
                    freq_data['sift'].append('')
                    freq_data['OMIM'].append('No entries')
                    freq_data['codons'].append('')
                    freq_data['codonsFraction'].append('')
                    freq_data['codonsUsageBias'].append('')
                    freq_data['gnomAD_gene_pLi'].append('')
                    freq_data['gnomAD_gene_pRec'].append('')
                    freq_data['gnomAD_gene_pNull'].append('')
                    freq_data['gnomAD_gene_synZ'].append('')
                    freq_data['gnomAD_gene_misZ'].append('')
                    freq_data['gnomAD_gene_loeuf'].append('')
                if 'phylopScore' in variant_dict:
                    freq_data['phyloP'].append(variant_dict['phylopScore'])
                else:
                    freq_data['phyloP'].append('')
                if 'regBase_prediction_v1' in variant_dict:
                    freq_data['RegBaseV1_RegPhred'].append(variant_dict['regBase_prediction_v1']['REG_PHRED'])
                    freq_data['RegBaseV1_PatPhred'].append(variant_dict['regBase_prediction_v1']['PAT_PHRED'])
                else:
                    freq_data['RegBaseV1_RegPhred'].append('')
                    freq_data['RegBaseV1_PatPhred'].append('')
                if 'CADDv1.6' in variant_dict:
                    freq_data['CADD_raw'].append(variant_dict['CADDv1.6']['CADD_raw'])
                    freq_data['CADD_PHRED'].append(variant_dict['CADDv1.6']['CADD_PHRED'])
                else:
                    freq_data['CADD_raw'].append('')
                    freq_data['CADD_PHRED'].append('')
                if 'regBase_v1.1.1' in variant_dict:
                    freq_data['DANN'].append(variant_dict['regBase_v1.1.1']['DANN'])
                    if 'LINSIGHT' in variant_dict['regBase_v1.1.1']:
                        freq_data['LINSIGHT'].append(variant_dict['regBase_v1.1.1']['LINSIGHT'])
                    else:
                        freq_data['LINSIGHT'].append('')
                    if 'FATHMM-XF' in variant_dict['regBase_v1.1.1']:
                        freq_data['FATHMM-XF'].append(variant_dict['regBase_v1.1.1']['FATHMM-XF'])
                    else:
                        freq_data['FATHMM-XF'].append('')
                else:
                    freq_data['DANN'].append('')
                    freq_data['LINSIGHT'].append('')
                    freq_data['FATHMM-XF'].append('')
                if 'ncER_v2' in variant_dict:
                    freq_data['ncER_perc'].append(variant_dict['ncER_v2'][0]['ncER_perc'])
                else:
                    freq_data['ncER_perc'].append('')
                if 'NCBoost_v20190108' in variant_dict:
                    freq_data['NCBoost_score'].append(variant_dict['NCBoost_v20190108']['NCBoost_score'])
                    freq_data['NCBoost_chr_rank_perc'].append(variant_dict['NCBoost_v20190108']['NCBoost_chr_rank_perc'])
                else:
                    freq_data['NCBoost_score'].append('')
                    freq_data['NCBoost_chr_rank_perc'].append('')
                if 'spliceAI' in variant_dict:
                    freq_data['spliceAI_accGainScore'].append(spliceai[0])
                    freq_data['spliceAI_accGainDist'].append(spliceai[1])
                    freq_data['spliceAI_accLossScore'].append(spliceai[2])
                    freq_data['spliceAI_accLossDist'].append(spliceai[3])
                    freq_data['spliceAI_donGainScore'].append(spliceai[4])
                    freq_data['spliceAI_donGainDist'].append(spliceai[5])
                    freq_data['spliceAI_donLossScore'].append(spliceai[6])
                    freq_data['spliceAI_donLossDist'].append(spliceai[7])
                else:
                    freq_data['spliceAI_accGainScore'].append('')
                    freq_data['spliceAI_accGainDist'].append('')
                    freq_data['spliceAI_accLossScore'].append('')
                    freq_data['spliceAI_accLossDist'].append('')
                    freq_data['spliceAI_donGainScore'].append('')
                    freq_data['spliceAI_donGainDist'].append('')
                    freq_data['spliceAI_donLossScore'].append('')
                    freq_data['spliceAI_donLossDist'].append('')
                if 'revel' in variant_dict:
                    freq_data['revel'].append(variant_dict['revel']['score'])
                else:
                    freq_data['revel'].append('')
                if 'primateAI' in variant_dict:
                    freq_data['primateAI'].append(variant_dict['primateAI'][0]['scorePercentile'])
                else:
                    freq_data['primateAI'].append('')
                if 'clinvar' in variant_dict:
                    freq_data['clinvar_[Phen - Signif]'].append(', '.join([str(item) for item in clinvar]))
                else:
                    freq_data['clinvar_[Phen - Signif]'].append('')
                if 'regulatoryRegions' in variant_dict:
                    freq_data['regulatoryRegionID'].append(variant_dict['regulatoryRegions'][0]['id'])
                    freq_data['regulatoryRegionType'].append(variant_dict['regulatoryRegions'][0]['type'])
                else:
                    freq_data['regulatoryRegionID'].append('')
                    freq_data['regulatoryRegionType'].append('')
                if 'inLowComplexityRegion' in variant_dict:
                    freq_data['inLowComplexityRegion'].append(variant_dict['inLowComplexityRegion'])
                else:
                    freq_data['inLowComplexityRegion'].append('')
                if 'dbsnp' in variant_dict:
                    freq_data['dbSNP'].append(variant_dict['dbsnp'][0])
                else:
                    freq_data['dbSNP'].append('')
                if 'mitomap' in variant_dict:
                    freq_data['MITOMAP_diseases'].append(', '.join([str(item) for item in variant_dict['mitomap']['diseases']]))
                    freq_data['MITOMAP_clinSignificance'].append(', '.join([str(item) for item in variant_dict['mitomap']['clinicalSignificance']]))
                else:
                    freq_data['MITOMAP_diseases'].append('')
                    freq_data['MITOMAP_clinSignificance'].append('')



    freq_df = pd.DataFrame({key:pd.Series(value) for key, value in freq_data.items()})
    freq_df.to_csv('{}/{}_Nirvana_combined_rare_variants.tsv'.format(sample, sample), mode='a', sep='\t', index=False, header=False)
    for el in in_list:
        freq_df.iloc[el:el+1].to_csv('{}/{}_Nirvana_combined_rare_variants_genes_list.tsv'.format(folder, sample, sample), mode='a', sep='\t', index=False, header=False)
    del freq_data
    del freq_df
    print('Finished processing {} variants on {} in {}'.format(len(positions_l), chrom, datetime.now() - begin_time))


print(' ')
end = datetime.now()
print(end)
print('JSON files decoding and filtering command for sample {} completed.'.format(sample))
print('The analysis completed in {}.'.format(end-start))
sys.stdout.close()
