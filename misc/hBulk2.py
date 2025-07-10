## Generate Tables from PCA

from sklearn.decomposition import PCA, SparsePCA
import time
from scipy.sparse import triu, coo_matrix
from scipy import signal,stats
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import rankdata

celltype = 'B-NK'
chr = 1
data2 = np.load(f'hBulk/{celltype}_100k.npy')
print(data2.shape)
data3 = data2
print('dropout: ', np.size(data3[data3 != 0]) / data3.shape[1] ** 2)

#data3 = CG_res
ngene = data3.shape[0]
E1 = coo_matrix(data3)

ave = np.array([np.mean(np.diag(data3,k=i)) for i in range(ngene)])
#print(ave)
ave[ave == 0] = 1


E1.data = E1.data/ave[E1.col - E1.row]
E1 = E1.toarray()
E1 = np.triu(E1)
E1 = E1 + (np.triu(E1,k=1)).T

F2 = E1
row, col = np.diag_indices_from(F2)
#F2[row, col] = F2.diagonal(0) / 2
F2[row, col] = 0


CG_res = F2

#
C2 = np.corrcoef(CG_res)
C2[np.isnan(C2)]=0

Data = C2


pca = PCA(n_components=5, random_state=0, svd_solver='arpack')
comp = pca.fit_transform(Data)
comp1 = comp[:,0]
comp2 = comp[:,1]
comp3 = comp[:,2]
print('var: ', np.var(Data), np.mean(Data))
#print('diag:', np.diag(C2[142:,142:]))
print('explianed var: ', pca.explained_variance_ratio_)
print('sing. val: ', pca.singular_values_)
print(comp1.size)
n_components=3
# resolution = 100000
# #comp = -comp
# pca_df = pd.DataFrame(comp[:,:3])
# pca_df.columns = ["PC{}".format(i) for i in range(1, n_components + 1)]
# pca_df["chrom"] = 'chr7'
# pca_df["start"] = np.arange(0, len(pca_df) * resolution, resolution)
# pca_df["end"] = pca_df["start"] + resolution
# pca_df = pca_df[["chrom", "start", "end"] + ["PC{}".format(i) for i in range(1, n_components + 1)]]
#
# pca_df.to_csv('Blood_table.tsv',sep="\t",header=True,index=None)


#################################################
E1 = comp[:,0]

ref = np.load('hBulk/CpG_100k.npy')

if np.corrcoef(E1,ref)[0,1]<0:
    E1 = -E1

print(np.corrcoef(E1,ref)[0,1])
################################################

window = 5
E2 = signal.correlate(E1,np.ones(window)/window,mode='same')

print(np.corrcoef(E2,ref)[0,1])

Comp = [E1, E2]
K = 2
fig, axes = plt.subplots(nrows=K, ncols=1, gridspec_kw={'height_ratios': [0.5]*K},
                         figsize=(10, 2.25), sharey='row')

ngene = comp1.size
ind_ls = [628,670]
for i in range(K):
    ax = axes[i]
    sns.despine(bottom=True, left=True, ax=ax)
    x, y = np.arange(ngene), Comp[i]
    ax.fill_between(x, y, 0, where=y >= 0, facecolor='C3', interpolate=True)
    ax.fill_between(x, y, 0, where=y <= 0, facecolor='C0', interpolate=True)
    #ax.set_ylim([-value, value])
    for ind in ind_ls:
        ax.axvline(x=ind-30, color='b', linestyle='-')
    #ax.set_xlabel('PCA', loc="right")
    ax.set_xticks([])
    ax.set_yticks([])

plt.savefig(f'plot/Bulk_hCD_{celltype}_100k.pdf')
plt.show()

np.save(f'hBulk/Bulk_{celltype}_{chr}_PCA_corr_100k.npy', E1)
#####################################################
def Size_comp(vec):
    vec1 = vec.copy()
    vec1[vec1 > 0] = 1
    vec1[vec1 < 0] = 0
    vec2 = vec1[1:] - vec1[:-1]
    res = 1
    res_ls = []
    for s in vec2:
        if s == 0:
            res += 1
        else:
            res_ls.append(res)
            res = 1
    res_ls.append(res)
    return res_ls


def CP_num_comp(vec):
    vec1 = vec.copy()
    vec1[vec1 > 0] = 1
    vec1[vec1 < 0] = 0
    vec2 = vec1[1:] - vec1[:-1]
    return np.sum(np.abs(vec2))


E1 = comp[:,0]

print(np.median(Size_comp(E1)))
print(CP_num_comp(E1))

print(np.median(Size_comp(E2)))
print(CP_num_comp(E2))