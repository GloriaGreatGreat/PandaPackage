import pandas as pd
import numpy as np
# pd.options.mode.chained_assignment = None

datalist = [] # datalist is for storing all dataframes

for i in range (1, 4): # load all datas from three .txt files into datalist, classified by Chromsome types
    df = pd.read_csv('sample'+str(i)+'_depths.txt',
    sep="\t", names=['Chromosome', 'Position','Depth'], header=None, low_memory=False)
    df = df.drop(df.index[[0]])
    datalist.append(df.ix[df['Chromosome'] == 'chr1'])
    datalist.append(df.ix[df['Chromosome'] == 'chr2'])
    datalist.append(df.ix[df['Chromosome'] == 'chr3'])

'''
# Question 1: this is a different approach, for verification
unique = 0
for i in range (0, 3):
    c1 = pd.merge(datalist[i].loc[:,['Position']],datalist[i+3].loc[:,['Position']],on='Position',how='inner')
    c2 = pd.merge(datalist[i+6].loc[:,['Position']],datalist[i+3].loc[:,['Position']],on='Position',how='inner')
    c1['key'] = c1['Position'].astype(str)
    c2['key'] = c2['Position'].astype(str)
    temp = datalist[i+3]
    temp['key'] = temp['Position'].astype(str)
    m = temp[~temp.key.isin(c1.key)]
    m['key'] = m['Position'].astype(str)
    unique += len(m[~m.key.isin(c2.key)])
print unique
'''

'''
# Question 2: this is a different implementation, for verification
# common = 0
# for i in range (0, 3):
#    common += len(pd.merge((pd.merge(datalist[i].loc[:, ['Position']], datalist[i+3].loc[:, ['Position']],
#    on='Position', how='inner')), datalist[i+6].loc[:, ['Position']],on='Position', how='inner'))
# print(common)
'''

# Question 1 and 2
common = 0 # for storing common genomic positions among three files
c1 = 0 # for storing common genomic positions among Chromosome1 and Chromosome2
c2 = 0 # for storing common genomic positions among Chromosome2 and Chromosome3
for i in range (0, 3):
    common += len(pd.merge((pd.merge(datalist[i].loc[:, ['Position']], datalist[i+3].loc[:, ['Position']],
    on='Position', how='inner')), datalist[i+6].loc[:, ['Position']],on='Position', how='inner'))
    c1 += len(pd.merge(datalist[i].loc[:,['Position']],datalist[i+3].loc[:,['Position']],on='Position',how='inner'))
    c2 += len(pd.merge(datalist[i+6].loc[:,['Position']],datalist[i+3].loc[:,['Position']],on='Position',how='inner'))
# Using the idea that the (unique of chr2) = (total of chr2) -
# (common of chr1 and chr2) - (common of chr2 and chr3) + (common of chr1, chr2 and chr3)
print "Number of genomic positions that are unique to sample2_depth.txt: " + str(len(datalist[3]+datalist[4]+datalist[5])-c1-c2+common)
print "Number of genomic positions that in common between all three depth files: " + str(common)

# Question 3
total = []
for i in range(0, 3):
    count = 0
    for j in range(0, 3):
        for row in datalist[i+j*3].itertuples():
            count += float(row[3]) * 0.001
    total.append(count/3000)
print "The average read depth for chromosome 1 is " + str(total[0]) + ", for chromosome 2 is " + str(total[1]) + ", for chromosome 3 is " + str(total[2])

# Question 4
location = ''
index = -1
value = 0
for i in range (0, 3):
    frames = [datalist[i].loc[:, ['Position', 'Depth']],
    datalist[3+i].loc[:, ['Position', 'Depth']],
    datalist[6+i].loc[:, ['Position', 'Depth']]]
    result = pd.concat(frames)
    result['Depth'] = result['Depth'].astype('float')
    after_mean = result.groupby(['Position']).mean()
    current_max = after_mean['Depth'].max()
    if current_max > value:
        index = i+1
        value = current_max
        location = after_mean['Depth'].idxmax()
print "The genomic position with largest average depth is Chromosome " + str(index) + " and Position " + str(location)
