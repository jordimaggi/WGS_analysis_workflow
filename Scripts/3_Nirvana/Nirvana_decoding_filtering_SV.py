import pandas as pd
from json_loading_SV import *
from datetime import datetime
import os
import json
import gzip
import math
import sys
from openpyxl import Workbook
import itertools
from itertools import islice
import subprocess


version = '1.0'

#Creating a list of gene identifiers from the genes of interest list for eye disorders
genes_file = snakemake.input.genes_list
genes_df = pd.read_csv(genes_file, sep='\t')
genes_list = genes_df['gene symbol'].tolist()
genes_lyst = list(dict.fromkeys(genes_list))

sample = snakemake.wildcards.sample
sys.stdout = open(snakemake.log.log, 'w')
sys.stderr = sys.stdout

start = datetime.now()
print(start)
print("Starting annotated SV JSON file decoding and filtering command for sample {}".format(sample))
print(' ')

#Loading JSON data and creating chr indeces dictionaries and variants and genes data with the json_loading function
positions, genes, samples, IDs= json_loading(snakemake.input.sv, snakemake.input.genes_list, snakemake.wildcards.sample, snakemake.wildcards.folder)
print('Number of structural variant postions in VCF:', len(positions))



#Creating the omim and gnmoAD_scores dictionaries
omim_genes = {}
gnomAD_scores = {}

for i4 in range(len(genes)):
    gene_dict = json.loads(genes[i4])

    if 'gnomAD' in gene_dict:
        gnomAD_scores[gene_dict['name']] = gene_dict['gnomAD']

#Defining the variables that will hold the variants data before writing them to files
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', \
          'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

variants_field = 'variants'
gnomad_field = 'gnomad'

freq_data = {'inHouse_previous_comment': [], 'chromosome': [], 'position': [], 'refAllele': [], 'altAllele': [], 'svEnd': [], 'svLength': [], 'ID': [], \
         'breakendEventId': [], 'quality': [], 'splitReadCount': [], 'pairedEndReadCount': [], 'mappingQuality': [], \
         'totalDepth': [], 'genotype': [], 'genotypeQuality': [], \
         'variantType': [],  'gene(s)': [], 'transcript_consequence': [], 'hgsvc': [], 'hgsvp': [], 'hgvsg': [], 'TADA_and_ClassifyCNV_overlap': [], \
         'ClassifyCNV_ACMG': [], 'ClassifyCNV_Score': [], 'ClassifyCNV_Dosage-sensitive_genes': [], \
         'TADA_pathogenicityScore': [], 'TADA_pathogenicityLabel': [], 'TADA_affectedGenesCount': [], 'TADA_geneDistance': [], 'TADA_exonOverlap': [], \
         'TADA_affectedEnhancersCount': [], 'TADA_enhancerDistance': [], 'TADA_enhancerConservation': [], 'TADA_CTCFdistance': [], 'TADA_boundaryDistance': [], \
         'gnomAD_all': [], 'gnomAD_NFE': [], 'gnomAD_max': [], 'gnomAD_max_ethnicity': [],'gnomAD_controls': [], 'gnomAD_homCount': [], \
         'caller': [], 'clinvar_[Phen - Signif]': [], \
         'clinGen_entry': [], 'clinGen_Overlap': [], 'clinGen_variant': [], 'clinGen_Interpretation': [], 'clinGen_phenotypes': [], \
         'RegBaseV1_RegPhred': [], 'RegBaseV1_PatPhred': [], 'CADD_raw': [], 'CADD_PHRED': [], 'DANN': [], 'revel': [], 'primateAI': [], 'LINSIGHT': [], 'FATHMM-XF': [], 'ncER_perc': [], 'NCBoost_score': [], 'NCBoost_chr_rank_perc': [], \
         'regulatoryRegionID': [], 'regulatoryRegionType': [], 'inLowComplexityRegion': []}


begin_time = datetime.now()


#Checking if filtered table already exists; if so, old file needs to be deleted or renamed prior to repeating the procedure
if os.path.exists(snakemake.output.rare):
    os.remove(snakemake.output.rare)
    freq_df = pd.DataFrame({key:pd.Series(value) for key, value in freq_data.items()})
    freq_df.to_csv(snakemake.output.rare, sep='\t', index=False)
    exists = True

#If it does not exists, it creates it with the freq_data dictionary keys as column labels
else:
    #Creating file with header
    freq_df = pd.DataFrame({key:pd.Series(value) for key, value in freq_data.items()})
    freq_df.to_csv(snakemake.output.rare, sep='\t', index=False)
    exists = False


for i in range(len(positions)):
    if i % 10000 == 0:
        print('Processed {} variants in {}.'.format(i, datetime.now() - begin_time))
    position_dict = json.loads(positions[i])
    c = False
    if variants_field in position_dict:
        for variant_dict in position_dict[variants_field]:

            #Cleaning the 'transcripts' column
            if 'transcripts' in variant_dict:
                transcripts = [[],[],[],[],[], []]
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
                transcripts[0] = list(dict.fromkeys(transcripts[0]))


            #Cleaning the 'clinvar' column
            clinvar = []
            if 'clinvar' in variant_dict:
                for entry in variant_dict['clinvar']:
                    if 'phenotypes' not in entry:
                        clinvar.append('[No phenotype provided' + ' - ' + ', '.join([str(item) for item in entry['significance']]) + ']')
                    else:
                        clinvar.append('[' + ', '.join([str(item) for item in entry['significance']]) + ' - ' + ', '.join([str(item) for item in entry['significance']]) + ']')
            clinvar = list(dict.fromkeys(clinvar))


            sample_dict = position_dict['samples'][0]

            if gnomad_field in variant_dict:
                freq = variant_dict[gnomad_field]['allAf']
            elif 'oneKg' in variant_dict:
                freq = variant_dict['oneKg']['allAf']
            else:
                freq = 0
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

            #If gnomAD all below 1% or not present in gnomAD, add data to freq_data
            if freq < 0.01:
                empt = True
                overlap_SV = 0
                if 'InHouse_comment_fields_SV' in position_dict:
                    for a in range(len(position_dict['InHouse_comment_fields_SV'])):
                        if position_dict['position'] == position_dict['InHouse_comment_fields_SV'][a]['start']:
                            if 'end' in variant_dict and variant_dict['end'] == position_dict['InHouse_comment_fields_SV'][a]['end']:
                                freq_data['inHouse_previous_comment'].append(position_dict['InHouse_comment_fields_SV'][a]['inHouse_Comment'])
                                empt = False
                            elif 'reciprocalOverlap' in position_dict['InHouse_comment_fields_SV'][a]:
                                if position_dict['InHouse_comment_fields_SV'][a]['reciprocalOverlap'] > overlap_SV:
                                    overlap_SV = position_dict['InHouse_comment_fields_SV'][a]['reciprocalOverlap']
                                    max_i = a
                        elif 'reciprocalOverlap' in position_dict['InHouse_comment_fields_SV'][a]:
                            if position_dict['InHouse_comment_fields_SV'][a]['reciprocalOverlap'] > overlap_SV:
                                overlap_SV = position_dict['InHouse_comment_fields_SV'][a]['reciprocalOverlap']
                                max_i = a
                if empt == True and overlap_SV > 0.9:
                    freq_data['inHouse_previous_comment'].append('Overlap {}: {}'.format(overlap_SV, position_dict['InHouse_comment_fields_SV'][max_i]['inHouse_Comment']))
                    empt = False
                if empt == True:
                    freq_data['inHouse_previous_comment'].append('')
                freq_data['chromosome'].append(position_dict['chromosome'])
                freq_data['position'].append(int(position_dict['position']))
                freq_data['ID'].append(IDs[i])
                if 'svEnd' in position_dict:
                    freq_data['svEnd'].append(position_dict['svEnd'])
                else:
                    freq_data['svEnd'].append('')
                if 'svLength' in position_dict:
                    freq_data['svLength'].append(position_dict['svLength'])
                else:
                    if 'svEnd' in position_dict:
                        freq_data['svLength'].append(position_dict['svEnd']-position_dict['position'])
                    else:
                        freq_data['svLength'].append('')
                if 'breakendEventId' in position_dict:
                    freq_data['breakendEventId'].append(position_dict['breakendEventId'])
                else:
                    freq_data['breakendEventId'].append('')
                freq_data['refAllele'].append(position_dict['refAllele'])
                freq_data['altAllele'].append(', '.join([str(item) for item in position_dict['altAlleles']]))
                if 'quality' in position_dict:
                    freq_data['quality'].append(position_dict['quality'])
                else:
                    freq_data['quality'].append('')
                if 'splitReadCounts' in sample_dict:
                    freq_data['splitReadCount'].append(sample_dict['splitReadCounts'][0])
                else:
                    freq_data['splitReadCount'].append('')
                if 'pairedEndReadCounts' in sample_dict:
                    freq_data['pairedEndReadCount'].append(sample_dict['pairedEndReadCounts'][0])
                else:
                    freq_data['pairedEndReadCount'].append('')
                if 'totalDepth' in sample_dict:
                    freq_data['totalDepth'].append(sample_dict['totalDepth'])
                else:
                    freq_data['totalDepth'].append('')
                if 'genotype' in sample_dict:
                    freq_data['genotype'].append(sample_dict['genotype'])
                else:
                    freq_data['genotype'].append('')
                if 'genotypeQuality' in sample_dict:
                    freq_data['genotypeQuality'].append(sample_dict['genotypeQuality'])
                else:
                    freq_data['genotypeQuality'].append('')
                freq_data['variantType'].append(variant_dict['variantType'])
                if 'hgvsg' in variant_dict:
                    freq_data['hgvsg'].append(variant_dict['hgvsg'])
                else:
                    freq_data['hgvsg'].append('')
                if 'transcripts' in variant_dict:
                    freq_data['gene(s)'].append(', '.join([str(item) for item in transcripts[0]]))
                    transcripts[1] = list(dict.fromkeys(transcripts[1]))
                    freq_data['transcript_consequence'].append(', '.join([str(item) for item in transcripts[1]]))
                    freq_data['hgsvc'].append(', '.join([str(item) for item in transcripts[2]]))
                    freq_data['hgsvp'].append(', '.join([str(item) for item in transcripts[3]]))
                else:
                    freq_data['gene(s)'].append('')
                    freq_data['transcript_consequence'].append('')
                    freq_data['hgsvc'].append('')
                    freq_data['hgsvp'].append('')
                if 'TADA_and_ClassifyCNV_annotations' in position_dict:
                    max_overlap = 0
                    for e in range(len(position_dict['TADA_and_ClassifyCNV_annotations'])):
                        if 'reciprocalOverlap' in position_dict['TADA_and_ClassifyCNV_annotations'][e]:
                            overlap = position_dict['TADA_and_ClassifyCNV_annotations'][e]['reciprocalOverlap']
                        else:
                            continue
                        if overlap > max_overlap:
                            max_overlap = overlap
                            c = e
                    if c != False:
                        freq_data['TADA_and_ClassifyCNV_overlap'].append(max_overlap)
                        if 'ClassifyCNV_ACMG' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['ClassifyCNV_ACMG'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['ClassifyCNV_ACMG'])
                            freq_data['ClassifyCNV_Score'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['ClassifyCNV_Score'])
                        else:
                            freq_data['ClassifyCNV_ACMG'].append('')
                            freq_data['ClassifyCNV_Score'].append('')
                        if 'ClassifyCNV_Dosage-sensitive_genes' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['ClassifyCNV_Dosage-sensitive_genes'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['ClassifyCNV_Dosage-sensitive_genes'])
                        else:
                            freq_data['ClassifyCNV_Dosage-sensitive_genes'].append('')
                        if 'TADA_pathogenicityScore' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_pathogenicityScore'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_pathogenicityScore'])
                        else:
                            freq_data['TADA_pathogenicityScore'].append('')
                        if 'TADA_pathogenicityLabel' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_pathogenicityLabel'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_pathogenicityLabel'])
                        else:
                            freq_data['TADA_pathogenicityLabel'].append('')
                        if 'TADA_affectedGenesCount' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_affectedGenesCount'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_affectedGenesCount'])
                        else:
                            freq_data['TADA_affectedGenesCount'].append('')
                        if 'TADA_geneDistance' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_geneDistance'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_geneDistance'])
                        else:
                            freq_data['TADA_geneDistance'].append('')
                        if 'TADA_exonOverlap' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_exonOverlap'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_exonOverlap'])
                        else:
                            freq_data['TADA_exonOverlap'].append('')
                        if 'TADA_affectedEnhancersCount' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_affectedEnhancersCount'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_affectedEnhancersCount'])
                        else:
                            freq_data['TADA_affectedEnhancersCount'].append('')
                        if 'TADA_enhancerDistance' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_enhancerDistance'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_enhancerDistance'])
                        else:
                            freq_data['TADA_enhancerDistance'].append('')
                        if 'TADA_enhancerConservation' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_enhancerConservation'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_enhancerConservation'])
                        else:
                            freq_data['TADA_enhancerConservation'].append('')
                        if 'TADA_CTCFdistance' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_CTCFdistance'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_CTCFdistance'])
                        else:
                            freq_data['TADA_CTCFdistance'].append('')
                        if 'TADA_boundaryDistance' in position_dict['TADA_and_ClassifyCNV_annotations'][c]:
                            freq_data['TADA_boundaryDistance'].append(position_dict['TADA_and_ClassifyCNV_annotations'][c]['TADA_boundaryDistance'])
                        else:
                            freq_data['TADA_boundaryDistance'].append('')
                    else:
                        freq_data['ClassifyCNV_ACMG'].append('')
                        freq_data['ClassifyCNV_Score'].append('')
                        freq_data['ClassifyCNV_Dosage-sensitive_genes'].append('')
                        freq_data['TADA_pathogenicityScore'].append('')
                        freq_data['TADA_pathogenicityLabel'].append('')
                        freq_data['TADA_affectedGenesCount'].append('')
                        freq_data['TADA_geneDistance'].append('')
                        freq_data['TADA_exonOverlap'].append('')
                        freq_data['TADA_affectedEnhancersCount'].append('')
                        freq_data['TADA_enhancerDistance'].append('')
                        freq_data['TADA_enhancerConservation'].append('')
                        freq_data['TADA_CTCFdistance'].append('')
                        freq_data['TADA_boundaryDistance'].append('')
                else:
                    freq_data['ClassifyCNV_ACMG'].append('')
                    freq_data['ClassifyCNV_Score'].append('')
                    freq_data['ClassifyCNV_Dosage-sensitive_genes'].append('')
                    freq_data['TADA_pathogenicityScore'].append('')
                    freq_data['TADA_pathogenicityLabel'].append('')
                    freq_data['TADA_affectedGenesCount'].append('')
                    freq_data['TADA_geneDistance'].append('')
                    freq_data['TADA_exonOverlap'].append('')
                    freq_data['TADA_affectedEnhancersCount'].append('')
                    freq_data['TADA_enhancerDistance'].append('')
                    freq_data['TADA_enhancerConservation'].append('')
                    freq_data['TADA_CTCFdistance'].append('')
                    freq_data['TADA_boundaryDistance'].append('')
                if len(freq_data['caller']) < len(freq_data['position']):
                    freq_data['caller'].append('')
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
                if gnomad_field in variant_dict:
                    freq_data['gnomAD_all'].append(freq)
                    freq_data['gnomAD_max'].append(gnomad_max[0])
                    freq_data['gnomAD_max_ethnicity'].append(gnomad_max[1])
                    try:
                        freq_data['gnomAD_NFE'].append(gnomad_NFE)
                    except:
                        freq_data['gnomAD_NFE'].append('')
                    if 'controlsAllAf' in variant_dict[gnomad_field]:
                        freq_data['gnomAD_controls'].append(variant_dict[gnomad_field]['controlsAllAf'])
                    else:
                        freq_data['gnomAD_controls'].append('')
                    freq_data['gnomAD_homCount'].append(variant_dict[gnomad_field]['allHc'])
                else:
                    freq_data['gnomAD_all'].append('')
                    freq_data['gnomAD_max'].append(gnomad_max[0])
                    freq_data['gnomAD_max_ethnicity'].append('')
                    freq_data['gnomAD_NFE'].append('')
                    freq_data['gnomAD_homCount'].append('')
                    freq_data['gnomAD_controls'].append('')
                if 'clingen' in position_dict:
                    entry = []
                    overlap = []
                    variant = []
                    interpretation = []
                    phenotypes = []
                    for els in position_dict['clingen']:
                        if 'variantType' in els:
                            entry.append(els['variantType'])
                        if 'annotationOverlap' in els:
                            overlap.append(els['annotationOverlap'])
                        var = 'chr' + str(els['chromosome']) + ':' + str(els['begin']) + '_' + str(els['end'])
                        variant.append(var)
                        if 'clinicalInterpretation' in els:
                            interpretation.append(els['clinicalInterpretation'])
                        if 'phenotypes' in els:
                            phenotypes.append(els['phenotypes'][0])
                    entry = list(dict.fromkeys(entry))
                    overlap = list(dict.fromkeys(overlap))
                    variant = list(dict.fromkeys(variant))
                    phenotypes = list(dict.fromkeys(phenotypes))
                    freq_data['clinGen_entry'].append(', '.join([str(item) for item in entry]))
                    freq_data['clinGen_Overlap'].append(', '.join([str(item) for item in overlap]))
                    freq_data['clinGen_variant'].append(', '.join([str(item) for item in variant]))
                    freq_data['clinGen_Interpretation'].append(', '.join([str(item) for item in interpretation]))
                    freq_data['clinGen_phenotypes'].append(', '.join([str(item) for item in phenotypes]))
                else:
                    freq_data['clinGen_entry'].append('')
                    freq_data['clinGen_Overlap'].append('')
                    freq_data['clinGen_variant'].append('')
                    freq_data['clinGen_Interpretation'].append('')
                    freq_data['clinGen_phenotypes'].append('')




freq_df = pd.DataFrame({key:pd.Series(value) for key, value in freq_data.items()})
freq_df.to_csv(snakemake.output.rare, mode='a', sep='\t', index=False, header=False)
del freq_data
del freq_df
print('Finished processing {} variants in {}'.format(i, datetime.now() - begin_time))



print(' ')
end = datetime.now()
print(end)
print('JSON files decoding and filtering command for sample {} completed.'.format(sample))
print('The analysis completed in {}.'.format(end-start))
sys.stdout.close()
