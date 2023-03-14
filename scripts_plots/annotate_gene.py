import pandas as pd
import numpy as np


def annotate(chr_, start, end):
    result = []
    for i in range(0, len(chrom)):
        if str(chr_) == str(chrom[i]):
            if start <= float(s_[i]) and end >= float(e_[i]):
                result.append(name[i])
    return result, len(result)

df = pd.read_csv('real/AR/CRC01/inferred_prof.csv', sep = ',')
chr_ = pd.read_csv('real/CNV_data/SegFixed_withFormalName_usedInScience.xls', sep = '\t')
df['CHR'], df['start'], df['end'] = chr_['CHR'], chr_['START'], chr_['END']
df = df[df['node6'] != df['node9']]
df = df[['CHR','start','end', 'node6','node9']]


source = pd.read_csv("real/colorectal.csv")
name = source['Gene Symbol']
chrom = source['Chr']
s_ = source['start']
e_ = source['end']
chrom = [23 if x=='X' else x for x in chrom]
chrom = [24 if x=='Y' else x for x in chrom]
selected_chrom = [1,2,3,4,5,6,7,8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]


gene = []
end = 0
first = 0
for index, row in df.iterrows():
    chr_, s, e = int(row['CHR'][3:]), int(row['start']), int(row['end'])
    if first == 0:
        c_ = chr_
    if c_ in set(selected_chrom) or chr_ in set(selected_chrom):
        #print(chr_)
        if first == 0:
            #print('new events', num)
            start = s
            end = e
            #ad = num
        if first == 1:
            if s == end + 1:
                end = e
            else:
                #if num != ad and c_ == 9 and num == 3:http://localhost:8890/notebooks/Documents/research_tempruns/BD_output/real_gene.ipynb#
                print(c_, start, end)
                out, counter = annotate(c_, start, end)
                print(','.join(out))
                if counter !=0:
                    print('')
                for i in range(0, counter):
                    gene.append(out[i])
                start = s
                end = e
                c_ = chr_
                #ad = num

        first = 1
