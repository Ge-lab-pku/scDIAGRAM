## Compute Bulk 100kb

## Compute Bulk Mb

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

celltype_ls = ['Ex'+str(i) for i in range(1,17)]


meta_data = pd.read_csv('mBC_meta_data.csv',index_col=0)
chr = 2
print(meta_data.iloc[:5,:])
print(meta_data.shape)


ind = meta_data['celltype'].isin(celltype_ls)
meta_data1 = meta_data.loc[ind,:]

Cell_num = meta_data1.shape[0]

print('Cellnumber: ', Cell_num)
chr_data = pd.read_table('mm10_chr_size.txt',header=None)
resolution = 100000
length = chr_data.iloc[chr-1, 1]
L = length // resolution + 1
Res = np.zeros((L,L))
Res = Res[30:,30:]
for c in range(Cell_num):
    cell = meta_data1.iloc[c,:]
    #print(cell.name)
    cc = cell.name
    bb = cc.split('_')
    cell1 = bb[0]
    batch = cc[-6:]
    print(c, cell1, batch)

    lib1 = 'mBC_'+batch
    cell1 = cell1.replace(',', '_')
    file = '../' + lib1 + f'/chr{chr}/' + cell1 + '.pairs'

    Mat = np.zeros((L, L))
    with open(file, 'rt') as f:
        for line in f:
            bb = line.split("\t")
            source = int(bb[1]) // resolution
            end = int(bb[3]) // resolution
            if source >= L or end >= L:
                continue
            Mat[source, end] += 1
            Mat[end, source] += 1
    row, col = np.diag_indices_from(Mat)
    Mat[row, col] = 0
    data = Mat[30:, 30:]

    Res += data

#############################################
A = Res/Cell_num
np.save(f'HiC_mat/Bulk_Ex_100k_{chr}.npy', A)
