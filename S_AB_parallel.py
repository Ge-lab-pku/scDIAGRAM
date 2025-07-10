## Specify AB after MCMC
## parallel

from scipy import stats
import pandas as pd
from cooltools.lib import numutils
import argparse
import numpy as np
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy import signal
from scipy.stats import binom,rankdata,spearmanr
# from cooltools.lib import numutils
# from cooltools.eigdecomp import _phase_eigs
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Implement SpecifyAB algorithm')
parser.add_argument('--file', type=str)
parser.add_argument('--chr', type=int, default=1)
parser.add_argument('--MH_file', type=str)
parser.add_argument('--ncpus', type=int, default=40)
parser.add_argument('--cpg_file', type=str, default='reference/mm10.CpG.100000.txt')
parser.add_argument('--reference', type=str, default='reference/mm10_chr_size.txt')
parser.add_argument('--resolution', type=int, default=100000)
parser.add_argument('--contact_thre', type=float, default=99.9)
parser.add_argument('--output_prefix', type=str, default="Ncut")
parser.add_argument('--output_dir', type=str, default="DistAB")
parser.add_argument('--input_version', type=str, default='v1')
parser.add_argument('--species', type=str, default='mm')
parser.add_argument('--before', type=bool, default=True)
#parser.add_argument('--cpg_file', type=str, default='reference/mm10.CpG.100000.txt')
args = parser.parse_args()



chr = args.chr
MH_ls_tot = np.load(args.MH_file)
print(MH_ls_tot.shape[0])





# LL = meta_data1.shape[0]
# print('Cell number: ', LL)



def Hamming(E1, E2):
    score = 0
    a1 = (E1 >= 0)
    a2 = (E2 >= 0)
    score += np.sum(a1 * a2)
    a1 = (E1 <= 0)
    a2 = (E2 <= 0)
    score += np.sum(a1 * a2)
    return score / np.size(E1)
def CpG_annotate2(cpg_ls, m_ls):
    H1 = np.zeros(N)
    K = len(cpg_ls)
    source = 0
    ind_ls = m_ls + [N]
    for k in range(K):
        ind = ind_ls[k]
        H1[source:ind] = cpg_ls[k]
        source = ind
    return H1



MH_heat = []
MH_heat2 = []


cor_res = []
cor_res2 = []
RNA_res = []
CpG = pd.read_table(args.cpg_file, header=None)
chr1 = f'chr{chr}'
CpG1 = rankdata(CpG.iloc[:, 2]) / CpG.shape[0]
CpG.iloc[:, 2] = CpG1
CpG1 = CpG.loc[CpG[0] == chr1, :]
linear_cpg_vector = CpG1[2].to_numpy()
CpG = linear_cpg_vector-0.5

ref = CpG.copy()
N = ref.shape[0]

if args.species == 'hg':
    centro_pos = pd.read_table('reference/centro_pos.txt',header=None,index_col=0)
    before_pos = centro_pos.iloc[args.chr-1, 0]
    end_pos = centro_pos.iloc[args.chr-1, 1]


def comp_score(Mat, ind_ls):
    ind1 = ind_ls == 1
    S1 = np.sum(Mat[ind1, :][:, ind1])
    ind2 = ind_ls == 0
    S2 = np.sum(Mat[ind2, :][:, ind2])
    return S1 + S2


E1_heat = []
def normal_cut(Mat, ind_ls):
    ind1 = ind_ls == 1
    ind2 = ind_ls == 0
    S = np.sum(Mat[ind1, :][:, ind2])
    d = np.sum(Mat,axis=1)
    k1 = np.sum(d[ind1])
    k2 = np.sum(d[ind2])
    return S/k1/k2

#bulk_data = np.mean(Data_tot[:, 3:, 3:], axis=0)
# E1 = np.load(f'HiC_mat/Bulk_Ex_{chr}_PCA_corr_100k_v2.npy')
# # E1 = E1[30:]
# E1[np.isnan(E1)]=0

# bulk_data = np.load(f'HiC_mat/Bulk_Ex_100k_{chr}.npy')
# # bulk_data = bulk_data[30:,30:]
# N = bulk_data.shape[0]


chr_data = pd.read_table(args.reference,header=None)
resolution = args.resolution
length = chr_data.iloc[chr - 1, 1]
L = length // resolution + 1
Mat = np.zeros((L, L))

def Compute(index):
    print(index)

    if args.input_version == "v2":
        global L
        Mat = np.zeros((L, L))
        file = cell_ls[index]
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
    elif args.input_version == "v1":
        Mat = np.load(file)
        if len(Mat.shape)==2:
            row, col = np.diag_indices_from(Mat)
            Mat[row, col] = 0
        else:
            Mat = Mat[index,:,:]
            row, col = np.diag_indices_from(Mat)
            Mat[row, col] = 0
    else:
        print("Undefined input version, error...")
        return 0
    
    X = Mat

    if args.species == "hg":
        if args.before==True and before_pos!=0:
            X = X[:before_pos, :before_pos]
        elif args.before==True and before_pos==0:
            print('Centromere is on the end, no small arms before centromere.')
            return 0
        else:
            X = X[end_pos:, end_pos:]
    elif args.species == "mm":
        end_pos = 3000000//resolution
        X = X[end_pos:, end_pos:]

    data = X.copy()


    ind_ls = list(MH_ls_tot[index, :])

    MH_ls = ind_ls

    A2 = data
    N = data.shape[0]

    A = A2
    mask = A.sum(axis=0) > 0

    OE, _, _, _ = numutils.observed_over_expected(A, mask)

    contact_thre = args.contact_thre
    clip_percentile = contact_thre
    if np.quantile(OE, contact_thre / 100) == 0:
        OE[OE > 0] = 1
    else:
        OE = np.clip(OE, 0, np.percentile(OE, clip_percentile))

    OE[~mask, :] = 0
    OE[:, ~mask] = 0

    row, col = np.diag_indices_from(OE)
    OE[row, col] = 0

    Res2 = OE
    #Res2[Res2>0]=1
    source_i = 0
    score_tot = []
    ind_ls = ind_ls + [N]
    for ind_i in ind_ls:
        source_j = 0
        score_ls = []
        for ind_j in ind_ls:
            data_do = Res2[source_j:ind_j, source_i:ind_i]
            score = np.mean(data_do)
            score_ls.append(score)
            source_j = ind_j
        score_tot.append(score_ls)
        source_i = ind_i

    score_tot = np.array(score_tot)
    score_tot = np.clip(score_tot, 0, np.percentile(score_tot, 99.9))
    D = np.sum(score_tot,axis=1)
    D1 = np.diag(1/D**(0.5))
    D1[np.isinf(D1)] = 0
    score_tot2 = np.matmul(D1,np.matmul(score_tot, D1))
    score_tot2[np.isnan(score_tot2)] = 0

    eigenvalue, eigenvector = np.linalg.eig(score_tot2)
    eigenvalue1 = eigenvalue[np.argsort(eigenvalue)]
    eigenvector1 = eigenvector[:, np.argsort(eigenvalue)]


    eigval_f = eigenvalue1[eigenvalue1 < 0.98]
    eigvec_f = eigenvector1[:, eigenvalue1 < 0.98]
    print(eigval_f[-1])


    res = eigvec_f[:, -1].copy()



    res1 = res.copy()
    res1.sort()
    SSres = -np.inf
    for s in res1[1:-1]:
        sign = 1
        res2 = res.copy()
        res2[res2 >= s] = 1
        res2[res2 < s] = 0
        s_vec0 = CpG_annotate2(res2, list(ind_ls[:-1]))
        #T0 = Hamming(s_vec0 - 0.5, CpG)
        T0 = np.corrcoef(s_vec0, CpG)[0,1]
        #print(s, T)
        if T0 < 0:
            T0 = -T0
            res2 = 1-res2
            sign = -sign
        if T0 > SSres:
            SSres = T0
            ssres = s
            ssres2 = res2
            ss_sign = sign
    # print(ssres, np.min(np.abs(res)),eigenvalue[1])
    print(ssres, np.min(np.abs(res)), SSres)
    res = ssres2
    res1 = (eigenvector1[:, -2] - ssres + 0.005)*ss_sign



    s_vec = CpG_annotate2(res, list(ind_ls[:-1]))
    s_vec2 = CpG_annotate2(res1, list(ind_ls[:-1]))

    s1 = (1-eigval_f[-1])/np.sum(score_tot)

    return s_vec, s_vec2, eigval_f[-1], s1


if args.input_version == "v2":
    cell_file = args.file

    cell_ls = []
    with open(cell_file, 'rt') as f:
        count = 0
        for line in f:
            cell_ls.append(line[:-1])

    LL = len(cell_ls)
else:
    Mat = np.load(file)
    LL = Mat.shape[0]
    
with Pool(processes=args.ncpus) as pool:
    result = pool.map(Compute, range(LL))
print("Program finished!")

res_ls = [s[0] for s in result]
res_ls2 = [s[1] for s in result]
res_ls3 = [s[2] for s in result]
res_ls4 = [s[3] for s in result]

res_ls = np.array(res_ls)
np.save(f'{args.output_dir}/{chr1}/{args.output_prefix}{chr1}.npy', res_ls)

res_ls2 = np.array(res_ls2)
print(res_ls2.shape)
np.save(f'{args.output_dir}/{chr1}/{args.output_prefix}{chr1}_real.npy', res_ls2)

res_ls3 = np.array(res_ls3)
print(res_ls3.shape)
np.save(f'{args.output_dir}/{chr1}/{args.output_prefix}{chr1}_eig.npy', res_ls3)

res_ls4 = np.array(res_ls4)
print(res_ls4.shape)
np.save(f'{args.output_dir}/{chr1}/{args.output_prefix}{chr1}_Compsc.npy', res_ls4)