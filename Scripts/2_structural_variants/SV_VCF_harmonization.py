import pandas as pd
import datetime
import os
import sys
import subprocess


version = '1.0'

sample = snakemake.wildcards.sample
folder = snakemake.wildcards.folder
sys.stdout = open("{}/{}/Console_outputs/SV_VCF_harmonization.{}.log".format(folder, sample, sample), 'w')
sys.stderr = sys.stdout
currentDate = datetime.date.today()
date = currentDate.strftime("%Y%m%d")

start = datetime.datetime.now()
print(start)
print("Starting annotated SV VCF files harmonization and filtering command for sample {}".format(sample))
print(' ')

#Loading each VCF one by one into a dataframe, edit data, filter, write to new, edited VCF
server = '{}/{}/SV_calling'.format(folder, sample)

#Defining the fields to be written
format_fields = 'GT:GQ:PR:SR:DP'
samples_fields = ['GT', 'GQ', 'PR', 'SR', 'DP']
info_fields = ['IMPRECISE', \
              'SVTYPE', \
              'SVLEN', \
              'END', \
              'CIPOS', \
              'CIEND', \
              'MATEID', \
              'EVENT', \
              'HOMLEN', \
              'HOMSEQ', \
              'INSLEN', \
              'INSSEQ', \
              'LEFT_INSSEQ', \
              'RIGHT_INSSEQ']

with open('{}/Scripts/2_structural_variants/Decoding/SV_VCFs_header'.format(folder)) as f:
    common_header = f.read()

#Load Manta file
file = snakemake.input.manta_in
f = open(file, "r")
Lines = f.readlines()
for i in range(len(Lines)):
    line = Lines[i]
    if '#CHROM' in line:
        head = i
manta = pd.read_csv(file, sep='\t', header=head)

manta = manta[manta['FILTER']=='PASS'].reset_index()
manta['ID'] = "Manta_" + manta['ID'].str.split('Manta').str[1]
manta['FORMAT'] = format_fields
for i in range(manta.shape[0]):
    infos = {}
    formats = {}
    info_data = manta.iloc[i].at['INFO'].split(';')
    dp = '.'
    for el in info_data:
        if el == 'IMPRECISE':
            infos['IMPRECISE'] = 'IMPRECISE'
            continue
        if el.split('=')[0] == 'BND_DEPTH':
            dp = el.split('=')[1]
            continue
        if 'SVINS' in el:
            if '_' in el:
                field = el.split('=')[0].split('SV')[0] + el.split('=')[0].split('SV')[1]
                infos[field] = el.split('=')[1]
            else:
                infos[el.split('=')[0].split('SV')[1]] = el.split('=')[1]
            continue
        if el.split('=')[0] in info_fields:
            if el.split('=')[0] == 'MATEID':
                infos[el.split('=')[0]] = 'Manta_' + el.split('=')[1].split('Manta')[1]
            else:
                infos[el.split('=')[0]] = el.split('=')[1]
    INFO = ''
    for el in info_fields:
        if el in infos:
            if len(INFO) > 1:
                INFO += ';'
            if el == 'IMPRECISE':
                INFO += el
            else:
                INFO += el + '=' + infos[el]

    format_data = manta.iloc[i].at[sample].split(':')
    del format_data[3]
    del format_data[1]
    formats['GT'] = format_data[0]
    formats['GQ'] = format_data[1]
    formats['PR'] = format_data[2].split(',')[1]
    if len(format_data) == 4:
        formats['SR'] = format_data[3].split(',')[1]
    else:
        formats['SR'] = '.'
    formats['DP'] = dp
    FORMAT = ''
    for el in samples_fields:
        if FORMAT != '':
            FORMAT += ':' + formats[el]
        else:
            FORMAT = formats[el]
    manta.at[i,'INFO']=INFO
    manta.at[i,sample]=FORMAT

manta.drop(labels='index', axis=1, inplace=True)
header = """##fileformat=VCFv4.3
##fileDate={}
##source=Manta_GenerateSVCandidates_v1.6.0
""".format(date)
header += common_header

output_VCF = snakemake.output.manta_out
with open(output_VCF, 'w') as vcf:
    vcf.write(header)
manta.to_csv(output_VCF, sep="\t", mode='a', index=False)


#WhamG VCF
file = snakemake.input.whamg_in
f = open(file, "r")
Lines = f.readlines()
for i in range(len(Lines)):
    line = Lines[i]
    if '#CHROM' in line:
        head = i
wham = pd.read_csv(file, sep='\t', header=head)

wham = wham[wham['FILTER']=='PASS'].reset_index()
wham['A'] = wham['INFO'].str.split(';').str[0].str.split('=').str[1].astype(int)
#Filter out variants with less than 5 supporting reads (as suggested on https://github.com/zeeev/wham)
wham = wham[wham['A']>5].reset_index()
wham.drop(['A', 'index', 'level_0'], axis=1, inplace=True)
wham['ID'] = "Wham"
wham['FORMAT'] = format_fields
for i in range(wham.shape[0]):
    infos = {}
    formats = {}
    info_data = wham.iloc[i].at['INFO'].split(';')
    dp = '.'
    for el in info_data:
        if el.split('=')[0] == 'SR':
            sr = el.split('=')[1]
            continue
        if el.split('=')[0] == 'TF':
            pe = el.split('=')[1]
            continue
        if el.split('=')[0] in info_fields:
            infos[el.split('=')[0]] = el.split('=')[1]
    INFO = ''
    for el in info_fields:
        if el in infos:
            if len(INFO) > 1:
                INFO += ';' + el + '=' + infos[el]
            else:
                INFO += el + '=' + infos[el]

    format_data = wham.iloc[i].at[sample].split(':')
    #for i2 in range(len(format_data)):
    #    formats[samples_fields[i2]] = format_data[i2]
    formats['GT'] = '.'
    formats['SR'] = sr
    formats['PR'] = pe
    formats['GQ'] = '.'
    formats['DP'] = '.'
    FORMAT = ''
    for el in samples_fields:
        if FORMAT != '':
            FORMAT += ':' + formats[el]
        else:
            FORMAT = formats[el]
    wham.at[i,'INFO']=INFO
    wham.at[i,sample]=FORMAT

header = """##fileformat=VCFv4.3
##fileDate={}
##source=WHAM-GRAPHENING_v1.7.0-311-g4e8c-dirty
""".format(date)
header += common_header

output_VCF = snakemake.output.whamg_out
with open(output_VCF, 'w') as vcf:
    vcf.write(header)
wham.to_csv(output_VCF, sep="\t", mode='a', index=False)


#Lumpy VCF
file = snakemake.input.lumpy_in
f = open(file, "r")
Lines = f.readlines()
for i in range(len(Lines)):
    line = Lines[i]
    if '#CHROM' in line:
        head = i
lumpy = pd.read_csv(file, sep='\t', header=head)

lumpy = lumpy[lumpy['QUAL']>0].reset_index()
lumpy = lumpy[lumpy['#CHROM']!='chrUn_gl000237'].reset_index()
lumpy['FILTER'] = 'PASS'
lumpy.drop(labels=['index', 'level_0'], axis=1, inplace=True)
lumpy = lumpy[lumpy['#CHROM']!='chr4_gl383529_alt'].reset_index()
lumpy.drop(labels=['index'], axis=1, inplace=True)
lumpy['ID'] = "Lumpy_" + lumpy['ID'].str[:]
lumpy['FORMAT'] = format_fields
for i in range(lumpy.shape[0]):
    infos = {}
    formats = {}
    info_data = lumpy.iloc[i].at['INFO'].split(';')
    dp = '.'
    for el in info_data:
        if el == 'IMPRECISE':
            infos['IMPRECISE'] = 'IMPRECISE'
            continue
        if el.split('=')[0] in info_fields:
            if el.split('=')[0] == 'MATEID':
                infos[el.split('=')[0]] = 'Lumpy_' + el.split('=')[1]
            else:
                infos[el.split('=')[0]] = el.split('=')[1]
    INFO = ''
    for el in info_fields:
        if el in infos:
            if len(INFO) > 1:
                INFO += ';'
            if el == 'IMPRECISE':
                INFO += el
            else:
                INFO += el + '=' + infos[el]

    format_data = lumpy.iloc[i].at[sample].split(':')
    formats['GT'] = format_data[0]
    formats['GQ'] = format_data[4]
    formats['PR'] = format_data[1]
    formats['SR'] = format_data[2]
    formats['DP'] = format_data[7]
    FORMAT = ''
    for el in samples_fields:
        if FORMAT != '':
            FORMAT += ':' + formats[el]
        else:
            FORMAT = formats[el]
    lumpy.at[i,'INFO']=INFO
    lumpy.at[i,sample]=FORMAT

header = """##fileformat=VCFv4.3
##fileDate={}
##source=LUMPY_v0.2.13
""".format(date)
header += common_header

output_VCF = snakemake.output.lumpy_out
with open(output_VCF, 'w') as vcf:
    vcf.write(header)
lumpy.to_csv(output_VCF, sep="\t", mode='a', index=False)



#Delly_sv VCF
file = snakemake.input.delly_sv_in
f = open(file, "r")
Lines = f.readlines()
for i in range(len(Lines)):
    line = Lines[i]
    if '#CHROM' in line:
        head = i
delly_sv = pd.read_csv(file, sep='\t', header=head)

delly_sv = delly_sv[delly_sv['FILTER']=='PASS'].reset_index()
delly_sv = delly_sv[delly_sv['#CHROM']!='chrUn_gl000237'].reset_index()
delly_sv.drop(labels=['index', 'level_0'], axis=1, inplace=True)
delly_sv = delly_sv[delly_sv['#CHROM']!='chr4_gl383529_alt'].reset_index()
delly_sv.drop(labels=['index'], axis=1, inplace=True)
delly_sv['ID'] = "Delly_" + delly_sv['ID']
delly_sv['FORMAT'] = format_fields
for i in range(delly_sv.shape[0]):
    infos = {}
    formats = {}
    info_data = delly_sv.iloc[i].at['INFO'].split(';')
    sr = '.'
    pe = '.'
    for el in info_data:
        if el.split('=')[0] == 'SR':
            sr = el.split('=')[1]
            continue
        if el.split('=')[0] == 'PE':
            pe = el.split('=')[1]
            continue
        if el == 'IMPRECISE':
            infos['IMPRECISE'] = 'IMPRECISE'
            continue
        if el.split('=')[0] in info_fields:
            if el.split('=')[0] == 'CHR2':
                infos['MATEID'] = 'Delly_' + el.split('=')[1]
            if el.split('=')[0] == 'POS2':
                infos['MATEID'] = infos['MATEID'] + '_' + el.split('=')[1]
            else:
                infos[el.split('=')[0]] = el.split('=')[1]
    if 'MATEID' in infos:
        infos['EVENT'] = infos['MATEID']
    INFO = ''
    for el in info_fields:
        if el in infos:
            if len(INFO) > 1:
                INFO += ';'
            if el == 'IMPRECISE':
                INFO += el
            else:
                INFO += el + '=' + infos[el]

    format_data = delly_sv.iloc[i].at[sample].split(':')
    formats['GT'] = format_data[0]
    formats['GQ'] = format_data[2]
    formats['PR'] = pe
    formats['SR'] = sr
    formats['DP'] = '.'
    FORMAT = ''
    for el in samples_fields:
        if FORMAT != '':
            FORMAT += ':' + formats[el]
        else:
            FORMAT = formats[el]
    delly_sv.at[i,'INFO']=INFO
    delly_sv.at[i,sample]=FORMAT

header = """##fileformat=VCFv4.3
##fileDate={}
##source=delly_v1.0.3
""".format(date)
header += common_header

output_VCF = snakemake.output.delly_sv_out
with open(output_VCF, 'w') as vcf:
    vcf.write(header)
delly_sv.to_csv(output_VCF, sep="\t", mode='a', index=False)


#Delly CNV VCF
file = snakemake.input.delly_cnv_in
f = open(file, "r")
Lines = f.readlines()
for i in range(len(Lines)):
    line = Lines[i]
    if '#CHROM' in line:
        head = i
delly_cnv = pd.read_csv(file, sep='\t', header=head)

delly_cnv = delly_cnv[delly_cnv['FILTER']=='PASS'].reset_index()
delly_cnv['ID'] = "Delly_" + delly_cnv['ID']
delly_cnv['FORMAT'] = format_fields
delly_cnv.drop(['index'], axis=1, inplace=True)
for i in range(delly_cnv.shape[0]):
    infos = {}
    formats = {}
    info_data = delly_cnv.iloc[i].at['INFO'].split(';')
    sr = '.'
    pe = '.'
    for el in info_data:
        if el == 'IMPRECISE':
            infos['IMPRECISE'] = 'IMPRECISE'
            continue
        if el.split('=')[0] in info_fields:
            infos[el.split('=')[0]] = el.split('=')[1]
    INFO = ''
    for el in info_fields:
        if el in infos:
            if len(INFO) > 1:
                INFO += ';'
            if el == 'IMPRECISE':
                INFO += el
            else:
                INFO += el + '=' + infos[el]

    format_data = delly_cnv.iloc[i].at[sample].split(':')
    gt = './.'
    if int(format_data[1]) < 2:
        delly_cnv.at[i,'ALT']='<DEL>'
        if int(format_data[1]) == 1:
            gt = '0/1'
        elif int(format_data[1]) == 0:
            gt = '1/1'
    elif int(format_data[1]) > 2:
        delly_cnv.at[i,'ALT']='<DUP>'
        if int(format_data[1]) == 3:
            gt = '0/1'
        elif int(format_data[1]) == 4:
            gt = '1/1'
    formats['GT'] = gt
    formats['GQ'] = format_data[3]
    formats['PR'] = pe
    formats['SR'] = sr
    formats['DP'] = '.'
    FORMAT = ''
    for el in samples_fields:
        if FORMAT != '':
            FORMAT += ':' + formats[el]
        else:
            FORMAT = formats[el]
    delly_cnv.at[i,'INFO']=INFO
    delly_cnv.at[i,sample]=FORMAT

header = """##fileformat=VCFv4.3
##fileDate={}
##source=delly_v1.0.3
""".format(date)
header += common_header

output_VCF = snakemake.output.delly_cnv_out
with open(output_VCF, 'w') as vcf:
    vcf.write(header)
delly_cnv.to_csv(output_VCF, sep="\t", mode='a', index=False)


#CNVpytor VCF
file = snakemake.input.cnvpytor_in
f = open(file, "r")
Lines = f.readlines()
for i in range(len(Lines)):
    line = Lines[i]
    if '#CHROM' in line:
        head = i
cnvpytor = pd.read_csv(file, sep='\t', header=head)

cnvpytor = cnvpytor[cnvpytor['FILTER']=='PASS'].reset_index()
cnvpytor['ID'] = "CNVpytor_" + cnvpytor['ID'].str.split('_').str[2] + '_' + cnvpytor['ID'].str.split('_').str[3]
cnvpytor['FORMAT'] = format_fields
cnvpytor.columns = [*cnvpytor.columns[:-1], sample]
cnvpytor.drop(['index'], axis=1, inplace=True)
for i in range(cnvpytor.shape[0]):
    infos = {}
    formats = {}
    info_data = cnvpytor.iloc[i].at['INFO'].split(';')
    sr = '.'
    pe = '.'
    for el in info_data:
        if el == 'IMPRECISE':
            infos['IMPRECISE'] = 'IMPRECISE'
            continue
        if el.split('=')[0] in info_fields:
            infos[el.split('=')[0]] = el.split('=')[1]
    INFO = ''
    for el in info_fields:
        if el in infos:
            if len(INFO) > 1:
                INFO += ';'
            if el == 'IMPRECISE':
                INFO += el
            else:
                INFO += el + '=' + infos[el]

    format_data = cnvpytor.iloc[i].at[sample].split(':')
    formats['GT'] = format_data[0]
    formats['GQ'] = '.'
    formats['PR'] = '.'
    formats['SR'] = '.'
    formats['DP'] = '.'
    FORMAT = ''
    for el in samples_fields:
        if FORMAT != '':
            FORMAT += ':' + formats[el]
        else:
            FORMAT = formats[el]
    cnvpytor.at[i,'INFO']=INFO
    cnvpytor.at[i,sample]=FORMAT


header = """##fileformat=VCFv4.3
##fileDate={}
##source=CNVpytor_v1.2.1
""".format(date)
header += common_header

output_VCF = snakemake.output.cnvpytor_out
with open(output_VCF, 'w') as vcf:
    vcf.write(header)
cnvpytor.to_csv(output_VCF, sep="\t", mode='a', index=False)



print(' ')
end = datetime.datetime.now()
print(end)
print('SV VCF files harmonization and filtering command for sample {} completed.'.format(sample))
print('The analysis completed in {}.'.format(end-start))
sys.stdout.close()
