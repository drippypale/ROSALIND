'''
Drippypale
drippypale@gmail.com

Problem:
    name: Edit Distance Alignment
    link: https://rosalind.info/problems/edta/

'''

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

seq1, seq2 = list(FASTA.read_seq_list('edta.in').values())
n, m = len(seq1), len(seq2)

import numpy as np

dp = np.zeros(shape=(n + 1, m + 1), dtype=np.int64)
dp[1:, 0], dp[0, 1:] = range(0, -n, -1), range(0, -m, -1)

bt = np.zeros(shape=(n + 1, m + 1), dtype=np.int64)
bt[1:, 0], bt[0, 1:] = 1, 2
bt_ = {
    0: (-1, -1), # match/mismatch
    1: (-1, 0), # del
    2: (0, -1)
}

for i in range(1, n + 1):
    for j in range(1, m + 1):
        if seq1[i - 1] == seq2[j - 1]:
            bt[i, j], dp[i, j] = 0, dp[i - 1][j - 1]
        else:
            bt[i, j], dp[i, j] = max((1, dp[i - 1, j] - 1), (2, dp[i, j - 1] - 1), (0, dp[i - 1, j - 1] - 1), key=lambda x: x[1])

print(abs(dp[n, m]))

seq1_, seq2_ = '', ''
i, j = n, m
i_, j_ = 0, 0
while(True):
    if bt[i, j] == 0:
        seq1_ += seq1[i - 1]
        seq2_ += seq2[j - 1]
    elif bt[i, j] == 1:
        seq2_ += '-'
        seq1_ += seq1[i - 1]
    else:
        seq1_ += '-'
        seq2_ += seq2[j - 1]

    if i == 1 and j == 1:
        break

    i_ = i + bt_[bt[i, j]][0]
    j_ = j + bt_[bt[i, j]][1]

    # print(i, j, bt_[bt[i, j]])

    i, j = i_, j_


print(seq1_[::-1])
print(seq2_[::-1])
        
    

