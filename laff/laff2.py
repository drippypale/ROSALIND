'''
Drippypale
drippypale@gmail.com

Problem:
    name: Local Alignment with Affine Gap Penalty
    link: https://rosalind.info/problems/laff/

'''


from utils.score_matrix import BLOSUM62

'''
I implemented a module to read the score matrix through the file,
returning the matrix alphabet index and the matrix itself as an Numpy ndArray
'''

BLOSUM_index, BLOSUM = BLOSUM62

a, b = -11, -1 # gap penalty: a + b(L - 1) (affine penalty)

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

p1, p2 = list(FASTA.read_seq_list('laff.in').values()) # Protein1 and Protein2 sequences
n, m = len(p1), len(p2)

import numpy as np
import math

'''
dp[0]: match/mismatch
dp[1]: gap in the first protein seq
dp[2]: gap in the second protein seq
'''
dp = np.zeros(shape=(3, n + 1, m + 1), dtype=np.int64)

# dp[1, 1:, 0], dp[1, 0, 1:] = range(a, n * b + a, b), 0
# dp[2, 0, 1:], dp[2, 1:, 0] = range(a, m * b + a, b), 0

bt = np.zeros(shape=(n + 1, m + 1), dtype=np.int64)

# bt[1, 1:, 0], bt[1, 0, 1:] = 11, 31
# bt[2, 0, 1:], bt[2, 1:, 0] = 22, 32

# bt[0, 0, 1:], bt[0, 1:, 0] = 3, 3

bt_ = {
    0: (-1, -1), # match/mismatch
    1: (-1, 0), # gap in the first protein seq
    2: (0, -1), # gap in the second protein seq
    3: (0, 0), # free ride from the start(0, 0) middle matrix

}

'''
dp and bt arrays will be saved to the npy files to load it quickly once they get generated for the first time.
this is useful for debuging the alignment part easily.
'''

try:
    with open('dp.npy', 'rb') as f1, open('bt.npy', 'rb') as f2:
        dp, bt = np.load(f1), np.load(f2)
except IOError:
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            dp[1, i, j] = max(
                dp[1, i - 1, j] + b, # extend the gap
                dp[0, i - 1, j] + a, # start the gap
            )
            dp[2, i, j] = max(
                dp[2, i, j - 1] + b, # extend the gap
                dp[0, i, j - 1] + a, # start the gap
            )
            bt[i, j], dp[0, i, j] = max(
                (3, 0), # free ride
                (0, dp[0, i - 1, j - 1] + BLOSUM[BLOSUM_index[p1[i - 1]], BLOSUM_index[p2[j - 1]]]),
                (1, dp[1, i, j]),
                (2, dp[2, i, j]),
                key=lambda x: x[1]
            )
    with open('dp.npy', 'wb') as f1, open('bt.npy', 'wb') as f2:
        np.save(f1, dp)
        np.save(f2, bt)

max_i, max_j = np.unravel_index(dp[0, :, :].argmax(), dp[0, :, :].shape)

print(dp[0, max_i, max_j])

i, j = max_i, max_j
i_, j_ = 0, 0
p1_, p2_ = '', ''

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