'''
Drippypale
drippypale@gmail.com

Problem:
    name: Local Alignment with Scoring Matrix
    link: https://rosalind.info/problems/loca/

'''


from utils.score_matrix import PAM250

'''
I implemented a module to read the score matrix through the file,
returning the matrix alphabet index and the matrix itself as an Numpy ndArray
'''

PAM_index, PAM = PAM250

g = -5 # gap penalty: gL (linear penalty)

from utils.helper import FASTA

'''
I implemented a helper module to do some common tasks like
reading a FASTA format input file, lcs, kmp, ...
which you can see the FASTA class below:

class FASTA:
    def read_seq_list(path: str) -> dict:
        seq_dict = dict()
        with open(path, 'r') as f:
            tag, seq = '', ''
            for l in f.read().splitlines():
                if l[0] == '>':
                    seq_dict[tag] = seq
                    tag, seq = l, ''
                else:
                    seq += l
            seq_dict[tag] = seq
        seq_dict.pop('')
        return seq_dict
'''

p1, p2 = list(FASTA.read_seq_list('loca.in').values()) # Protein1 and Protein2 sequences
n, m = len(p1), len(p2)

import numpy as np

dp = np.zeros(shape=(n + 1, m + 1), dtype=np.int64)
dp[1:, 0], dp[0, 1:] = range(0, -5 * n, -5), range(0, -5 * m, -5)

bt = np.zeros(shape=(n + 1, m + 1), dtype=np.int64)
bt[1:, 0], bt[0, 1:] = 1, 2
bt_ = {
    0: (-1, -1), # match/mismatch
    1: (-1, 0), # del
    2: (0, -1), # ins
    3: (0, 0) # free ride from the start (0, 0) !
}

for i in range(1, n + 1):
    for j in range(1, m + 1):
        if p1[i - 1] == p2[j - 1]:
            bt[i, j], dp[i, j] = 0, dp[i - 1, j - 1] + PAM[PAM_index[p1[i - 1]], PAM_index[p1[i - 1]]]
        else:
            bt[i, j], dp[i, j] = max(
                (3, dp[1, 1]),
                (1, dp[i - 1, j] + g),
                (2, dp[i, j - 1] + g),
                (0, dp[i - 1, j - 1] + PAM[PAM_index[p1[i - 1]], PAM_index[p2[j - 1]]]),
                key=lambda x: x[1]
            )

max_i, max_j = np.unravel_index(dp.argmax(), dp.shape)
print(dp[max_i, max_j])

p1_, p2_ = '', ''
i, j = max_i, max_j
i_, j_ = 0, 0
while True:
    if bt[i, j] == 0:
        p1_ += p1[i - 1]
        p2_ += p2[j - 1]
    elif bt[i, j] == 1:
        p1_ += p1[i - 1]
    elif bt[i, j] == 2:
        p2_ += p2[j - 1]
    else:
        break

    if i == 1 and j == 1:
        break

    i_ = i + bt_[bt[i, j]][0]
    j_ = j + bt_[bt[i, j]][1]

    i, j = i_, j_

print(p1_[::-1])
print(p2_[::-1])
