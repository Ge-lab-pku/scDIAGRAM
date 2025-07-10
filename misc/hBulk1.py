import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

celltype = 'B-NK'
meta_data = pd.read_csv('PseudoTime/diffMap2.csv',index_col=0)
print(meta_data.iloc[:5,:])
LL = meta_data.shape[0]
print('Cellnumber: ', LL)



cell_ls = list(meta_data.index)
resolution = 100000
chr = 1

chr = f'chr{chr}'
CpG = pd.read_table('reference/hg19.100kbin.CpG.txt', header=0)
CpG1 = CpG.loc[CpG['#1_usercol'] == chr, :]
CpG1['ratio'] = (CpG1['13_user_patt_count'] / (CpG1['12_seq_len'] - CpG1['10_num_N'])).fillna(0)
linear_cpg_vector = CpG1['ratio'].to_numpy()
ref = linear_cpg_vector.copy()

L = ref.shape[0]
np.save('hBulk/CpG_100k.npy',ref)
Res = np.zeros((L,L))

celltype_ls = list(meta_data['celltype'])
print(len(celltype_ls))
print(celltype_ls[:5])
count = 0
for index in range(LL):
    print(index)
    if celltype_ls[index]!=celltype:
        continue
    ss = cell_ls[index]
    bb = ss.split('_')
    lib = ss[-7:]
    cell = bb[0]
    cell1 = cell.replace(',', '_')
    file = '../hCD' + lib + f'/{chr}/' + cell1 + '.pairs'
    # global L


    Mat = np.zeros((L, L))
    with open(file, 'rt') as f:
        for line in f:
            bb = line.split("\t")
            source = int(bb[1]) // resolution
            target = int(bb[3]) // resolution
            if source >= L or target >= L:
                continue
            Mat[source, target] += 1
            Mat[target, source] += 1
    row, col = np.diag_indices_from(Mat)
    Mat[row, col] = 0


    Res += Mat
    count += 1

print(Res.shape)

Res = Res/count

val1 = np.quantile(Res,.9)
plt.matshow(Res,cmap='bwr',vmin=0,vmax=val1,aspect='auto')
plt.show()

np.save(f'hBulk/{celltype}_100k.npy',Res)


