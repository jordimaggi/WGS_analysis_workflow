
from datetime import datetime
import pandas as pd
import os
import glob

sample = snakemake.wildcards.sample
folder = snakemake.wildcards.folder
date = datetime.now().date()
genes_file = snakemake.input.genes_list
genes_file_name = genes_file.split('/')[-1][:-4]
sys.stdout = open(snakemake.log.log, 'w')

start = datetime.now()
print(start)
print("Starting merging of filtered tsv files for sample {} and genes in file {}".format(sample, genes_file))
print(' ')

if os.path.exists(snakemake.input.short):
    short_variants = pd.read_csv(snakemake.input.short, sep='\t')
    short_variants.insert(0, 'Comment', '')
else:
    print("File not found: {}".format(snakemake.input.short))

if os.path.exists(snakemake.input.sv):
    sv_variants = pd.read_csv(snakemake.input.sv, sep='\t')
    sv_variants.insert(0, 'Comment', '')
else:
    print("File not found: v".format(snakemake.input.sv))


variants = pd.concat([short_variants, sv_variants], ignore_index=True, sort=False)
variants.sort_values(by=['chromosome', 'position', 'refAllele', 'altAllele'], ignore_index=True, inplace=True)
cols = ['Comment', 'inHouse_previous_comment', 'reported', 'chromosome', 'position', 'ID', 'refAllele', 'altAllele', 'svEnd', 'svLength', 'breakendEventId', 'quality', 'splitReadCount', 'pairedEndReadCount', \
         'fisherStrandBias', 'mappingQuality', 'alleleDepths', 'variantFrequencies', 'totalDepth', 'genotype', 'genotypeQuality', 'caller', \
         'variantType', 'hgvsg', 'gene(s)', 'transcript_consequence', 'hgsvc', 'hgsvp', 'gnomAD_all', 'gnomAD_NFE', 'gnomAD_max', 'gnomAD_max_ethnicity', 'gnomAD_controls', 'gnomAD_homCount', \
         'clinvar_[Phen - Signif]', 'OMIM', 'codons', 'codonsFraction', 'codonsUsageBias', 'phyloP', \
         'RegBaseV1_RegPhred', 'RegBaseV1_PatPhred', 'CADD_raw', 'CADD_PHRED', 'DANN', 'revel', 'primateAI', 'LINSIGHT', 'FATHMM-XF', 'ncER_perc', 'NCBoost_score', 'NCBoost_chr_rank_perc', \
         'spliceAI_accGainScore', 'spliceAI_accGainDist', 'spliceAI_accLossScore', 'spliceAI_accLossDist', 'spliceAI_donGainScore', 'spliceAI_donGainDist', 'spliceAI_donLossScore', 'spliceAI_donLossDist', \
         'polyPhen', 'sift', 'TADA_and_ClassifyCNV_overlap', 'ClassifyCNV_ACMG', 'ClassifyCNV_Score', 'ClassifyCNV_Dosage-sensitive_genes', 'TADA_pathogenicityScore', 'TADA_pathogenicityLabel', \
         'TADA_affectedGenesCount', 'TADA_geneDistance', 'TADA_exonOverlap', 'TADA_affectedEnhancersCount', 'TADA_enhancerDistance', 'TADA_enhancerConservation', 'TADA_CTCFdistance', 'TADA_boundaryDistance', \
         'clinGen_entry', 'clinGen_Overlap', 'clinGen_variant', 'clinGen_Interpretation', 'clinGen_phenotypes', \
         'gnomAD_gene_pLi', 'gnomAD_gene_pRec', 'gnomAD_gene_pNull', 'gnomAD_gene_synZ', 'gnomAD_gene_misZ', 'gnomAD_gene_loeuf', \
         'regulatoryRegionID', 'regulatoryRegionType', 'inLowComplexityRegion', 'dbSNP', 'MITOMAP_diseases', 'MITOMAP_clinSignificance']

variants = variants[cols]
variants.to_csv(snakemake.output.rare, sep='\t', index=False)


if os.path.exists('{}/{}/Filtered_variants/{}_filtered_structural_and_short_variants.xlsx'.format(folder, sample, sample)):
    with pd.ExcelWriter('{}/{}/Filtered_variants/{}_filtered_structural_and_short_variants.xlsx'.format(folder, sample, sample), engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            variants.to_excel(writer, sheet_name='{}_variants'.format(genes_file_name), index=False)
else:
    if os.path.exists('{}/{}/Filtered_variants'.format(folder, sample)) == False:
        os.makedirs('{}/{}/Filtered_variants'.format(folder, sample))
    with pd.ExcelWriter('{}/{}/Filtered_variants/{}_filtered_structural_and_short_variants.xlsx'.format(folder, sample, sample), engine='openpyxl') as writer:
            variants.to_excel(writer, sheet_name='{}_variants'.format(genes_file_name), index=False)

print(' ')
end = datetime.now()
print(end)
print('Merging of filtered tsv files for sample {} for gene list {} completed.'.format(sample, genes_file_name))
print('The analysis completed in {}.'.format(end-start))
sys.stdout.close()
