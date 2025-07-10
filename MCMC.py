## Implement MCMC in the scDIAGRAM

import argparse
import time
from collections import Counter
from multiprocessing import Pool
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import signal, stats
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.preprocessing import quantile_transform
#from cooltools.lib import numutils

parser = argparse.ArgumentParser(description='Implement MH algorithm')

parser.add_argument('--file', type=str)
parser.add_argument('--index', type=int, default=0)
parser.add_argument('--chr', type=int, default=7)
parser.add_argument('--start', type=int)
parser.add_argument('--end', type=int)
parser.add_argument('--step', type=int, default=1)
parser.add_argument('--terminate', type=int, default=10000)
parser.add_argument('--out_file', type=str, default='single-cell_res.out')
#parser.add_argument('--CP_term', type=int, default=50)
parser.add_argument('--contact_thre', type=float, default=99.9)
#parser.add_argument('--cpg_file', type=str, default=None)
parser.add_argument('--before', type=bool, default=True)
parser.add_argument('--species', type=str, default='mm')
parser.add_argument('--resolution', type=int, default=100000)
parser.add_argument('--binary', type=int, default=False)
parser.add_argument('--input_version', type=str, default='v1')
parser.add_argument('--reference', type=str, default='reference/mm10_chr_size.txt')
parser.add_argument('--output', type=str, default='MCMC_res.npy')
#parser.add_argument('--out_file', type=str, default='single-cell_res.out')


args = parser.parse_args()
print(args)

file = args.file
start = args.start
end = args.end
step = args.step
index = args.index
outFile = args.out_file
#CP_term = args.CP_term

contact_thre = args.contact_thre
out_file = args.out_file


################################################################################
print('Preparation. Compute the variance ...')


chr = f'chr{args.chr}'
resolution = args.resolution

chr_data = pd.read_table('reference/mm10_chr_size.txt',header=None)
length = chr_data.iloc[args.chr - 1, 1]
L = length // resolution + 1

def cross_entropy(k, N):
    if k == 0 or k == N:
        return 0
    return k * np.log(k / N) + (N - k) * np.log(1 - k / N)



if args.species == 'hg':
    centro_pos = pd.read_table('reference/centro_pos.txt',header=None,index_col=0)
    before_pos = centro_pos.iloc[args.chr-1, 0]
    end_pos = centro_pos.iloc[args.chr-1, 1]


def comp(file, index=0):
    if args.input_version == "v2":
        global L
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
    # X = Mat[30:, 30:]



    X = Mat
    if args.binary == True:
        X[X > 0] = 1

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


    Z = X.copy()
    n = X.shape[0]
    print(n)

    C2 = X.copy()
    X = np.triu(X)
    Var = np.var(X[X != 0])
    print('var: ', Var)
    T = Var * 2

    time0 = time.time()
    logR = np.zeros((n, n))
    S_square = np.zeros((n + 1, n + 1))
    S = np.zeros((n + 1, n + 1))
    score1 = 0
    score2 = 0
    X_square = X ** 2
    # if T>0:
    for t in range(1, n + 1):
        time1 = time.time()
        # print(t,time1-time0)
        if t == 1:
            score1 += X[0, 0]
            score2 += X_square[0, 0]
            S[1, 1] = score1
            S_square[1, 1] = score2
            for s in range(2, n + 1):
                inc1 = np.sum(X[:s, s - 1])
                inc2 = np.sum(X_square[:s, s - 1])
                score1 += inc1
                score2 += inc2
                S[1, s] = score1
                S_square[1, s] = score2
        else:
            score1 = X[t - 2, t - 2]
            score2 = X_square[t - 2, t - 2]
            for s in range(t, n + 1):
                score1 += X[t - 2, s - 1]
                score2 += X_square[t - 2, s - 1]
                S[t, s] = S[t - 1, s] - score1
                S_square[t, s] = S_square[t - 1, s] - score2

            # X1 = X[(t - 1):s, (t - 1):s]

    time1 = time.time()
    print('used time: ', time1 - time0)

    N = X.shape[0]

    def logphi1(t, s):
        Nk = np.sum(X[t - 1:s, t - 1:s] != 0)
        N = (s - t + 1) * (s - t + 2) / 2
        if Nk == 0:
            return 0
        # N = (s - t + 1) * (s - t + 2) / 2
        V1 = S_square[t, s] / Nk - (S[t, s] / Nk) ** 2
        if T > 0:
            res = -Nk * V1 / T + cross_entropy(Nk, N)
        else:
            res = cross_entropy(S[t, s], N)
        return res

    def logphi2(t1, s1, t2, s2):
        Nk = np.sum(X[t1 - 1:s1, t2 - 1:s2] > 0)
        N = (s1 - t1 + 1) * (s2 - t2 + 1)
        if Nk == 0:
            return 0
        # N = (s1 - t1 + 1) * (s2 - t2 + 1)
        S1 = S[s1 + 1, t2 - 1] - S[t1, t2 - 1] - S[s1 + 1, s2] + S[t1, s2]
        S2 = S_square[s1 + 1, t2 - 1] - S_square[t1, t2 - 1] - S_square[s1 + 1, s2] + S_square[t1, s2]
        V1 = S2 / Nk - (S1 / Nk) ** 2
        if T > 0:
            res = -Nk * V1 / T + cross_entropy(Nk, N)
        else:
            res = cross_entropy(S1, N)
        return res

    def comp_score(m_ls, N):
        m_ls = [0] + m_ls + [N]
        kT = len(m_ls)
        score = 0
        for i in range(kT - 1):
            source = m_ls[i]
            end = m_ls[i + 1]
            score += logphi1(source + 1, end)

        for i in range(kT - 1):
            source1 = m_ls[i]
            end1 = m_ls[i + 1]
            for j in range(i + 1, kT - 1):
                source2 = m_ls[j]
                end2 = m_ls[j + 1]
                score += logphi2(source1 + 1, end1, source2 + 1, end2)

        return score




    Ind = np.arange(1, N)

    res_ls_tot = {}

    k_ls = np.arange(start, end, step)
    #print(k_ls)

    for k in k_ls:
        explode = 0

        count = 0

        np.random.seed(1234)
        m_ls = np.random.choice(N - 1, size=k) + 1
        m_ls.sort()
        m_ls = list(m_ls)
        ## T for variance*2.
        m_res = list(set(Ind) - set(m_ls))
        # print(m_ls)

        score = comp_score(m_ls, N)

        #L = {}

        logp = 0

        LT = []
        terminate = args.terminate
        res_ls = {}
        Likelih_ls = {}
        np.random.seed(1234)
        for t in range(args.end + 1):
            Likelih_ls[t] = -np.inf
        while True:
            m_ls_cand = m_ls.copy()
            m_ls_old = m_ls.copy()

            
            k = len(m_ls)
            flip1 = np.random.randint(N - k - 1)
            flip2 = np.random.randint(k)
            m_ls_cand.pop(flip2)
            m_ls_cand.append(m_res[flip1])
            m_ls_cand.sort()

            score2 = comp_score(m_ls_cand, N)

            prob_log = (score2 - score) + (len(m_ls_cand) - len(m_ls)) * logp

            accep_log = min(0, prob_log)
            trial_prob_log = np.log(np.random.rand())
            if trial_prob_log < accep_log:
                m_ls = m_ls_cand
                m_res = list(set(Ind) - set(m_ls))
                score = score2
                k = len(m_ls)
                if score > Likelih_ls[k]:
                    res_ls[k] = m_ls
                    Likelih_ls[k] = score

            if count % 100 == 0:
                print(count, m_ls, len(m_ls), accep_log, m_ls_cand, score)
            LT.append(len(m_ls))
            count += 1
            if count == terminate:
                break

        if explode == 1:
            print('explode...')
            continue


        k_most = k
        with open(outFile, 'a') as f:
            print(res_ls[k_most], Likelih_ls[k_most], k_most, file=f)
            print('No.', index, 'cell being processed on chr', args.chr, file=f)
        print(res_ls[k_most], Likelih_ls[k_most], k_most)

    #print(res_ls_tot)
    return res_ls[k_most]

res = comp(file)
print(res)

res = np.array(res)
np.save(args.output, res)
print(res.shape)