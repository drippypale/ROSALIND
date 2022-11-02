'''
Drippypale
drippypale@gmail.com

problem link: https://rosalind.info/problems/gaff/

'''

from multiprocessing.spawn import import_main_path
import numpy as np

import sys

BLOSUM_index: dict  # holds the index of each amino in the BLOSUM matrix
BLOSUM: np.ndarray  # BLOSUM matrix

with open('BLOSUM62', 'r') as f:
    lines = f.read().splitlines()
    BLOSUM_index = {c: i for i, c in enumerate(lines[0].split(' '))}
    n = len(BLOSUM_index)
    BLOSUM = np.zeros(shape=(n, n), dtype=np.int64)
    for i in range(1, n + 1):
        for j, j_ in enumerate(lines[i].split(' ')):
            BLOSUM[i - 1, j] = int(j_)

# print(BLOSUM_index)
# print(BLOSUM)

a, b = -11, -1 # gap penalty: a + b(L - 1)

from utils.helper import FASTA

p1, p2 = list(FASTA.read_seq_list('gaff.in').values()) # Protein1 and Protein2 sequences
n, m = len(p1), len(p2)

dp = np.zeros(shape=(n + 1, m + 1), dtype=np.int64)
dp[1:, 0], dp[0, 1:] = range(0, -n, -1), range(0, -m, -1)
dp[1:, 0] += -11
dp[0, 1:] += -11


bt = np.zeros(shape=(n + 1, m + 1), dtype=np.int64)
bt[1:, 0], bt[0, 1:] = 1, 2
bt_ = {
    0: (-1, -1), # match/mismatch
    1: (-1, 0), # del
    2: (0, -1)
}

for i in range(1, n + 1):
    for j in range(1, m + 1):
        if p1[i - 1] == p2[j - 1]:
            bt[i, j], dp[i, j] = 0, dp[i - 1, j - 1] + BLOSUM[BLOSUM_index[p1[i - 1]], BLOSUM_index[p1[i - 1]]]
        else:
            bt[i, j], dp[i, j] = max(
                (1, dp[i - 1, j] + (b if bt[i - 1, j] == 1 else a)),
                (2, dp[i, j - 1] + (b if bt[i, j - 1] == 2 else a)),
                (0, dp[i - 1, j - 1] + BLOSUM[BLOSUM_index[p1[i - 1]], BLOSUM_index[p2[j - 1]]]),
                key=lambda x: x[1]
            )
print(dp[n, m])

p1_, p2_ = '', ''
i, j = n, m
i_, j_ = 0, 0
while(True):
    if bt[i, j] == 0:
        p1_ += p1[i - 1]
        p2_ += p2[j - 1]
    elif bt[i, j] == 1:
        p2_ += '-'
        p1_ += p1[i - 1]
    else:
        p1_ += '-'
        p2_ += p2[j - 1]

    if i == 1 and j == 1:
        break

    i_ = i + bt_[bt[i, j]][0]
    j_ = j + bt_[bt[i, j]][1]

    i, j = i_, j_

print(p1_[::-1])
print(p2_[::-1])
