## Compute scCpG for a chr
## Test on chr1
## 100kb

import pandas as pd
import numpy as np
from scipy.stats import rankdata
from multiprocessing import Pool

CpG = pd.read_table('reference/hg19.100kbin.CpG.txt', header=0)
print(CpG.iloc[:5,:])
RRes = []

meta_data = pd.read_csv('PseudoTime/diffMap2.csv',index_col=0)

print(meta_data.shape)
print(meta_data.iloc[:5,:5])
LL = meta_data.shape[0]
resolution = 100000

c = 3

chr1 = f'chr{c}'
CpG1 = CpG.loc[CpG['#1_usercol'] == chr1, :]
CpG1['ratio'] = (CpG1['13_user_patt_count'] / (CpG1['12_seq_len'] - CpG1['10_num_N'])).fillna(0)
linear_cpg_vector = CpG1['ratio'].to_numpy()
ref = linear_cpg_vector.copy()

print('chr', c, ' cells: ', LL)
Res = []

N1 = ref.shape[0]
# N2 = Data.shape[1]-3
#
# if N1>N2:
#     ref = ref[:N2]

cell_ls = list(meta_data.index)
print(len(cell_ls))
print(cell_ls[:5])

# chr_data = pd.read_table('mm10_chr_size.txt',header=None)
L = ref.shape[0]
for cc in range(LL):
    print(c, cc)
    ss = cell_ls[cc]
    bb = ss.split('_')
    lib = ss[-7:]
    cell = bb[0]
    cell1 = cell.replace(',', '_')
    file = '../hCD' + lib + f'/chr{c}/' + cell1 + '.pairs'

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

    # data[data > 0] = 1
    # N2 = data.shape[0]
    #
    # N3 = min(N1, N2)
    # ref = ref[:N3]
    # data = data[:N3, :N3]

    s_vec = np.dot(Mat, ref) / np.sum(Mat, axis=0)
    s_vec[np.isnan(s_vec)] = ref[np.isnan(s_vec)]
    Res.append(s_vec)

Res = np.array(Res)
print(Res.shape)

np.save(f'hCpG1_{chr1}.npy', Res)
# RRes.append(Res)
