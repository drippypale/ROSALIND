'''
Drippypale
drippypale@gmail.com

Problem:
    name: Multiple Alignment
    link: https://rosalind.info/problems/mult/

'''

from utils.helper import FASTA
import numpy as np

seq_list = list(FASTA.read_seq_list('mult.in').values())

seq_list = ['$' + x for x in seq_list]
 
shape = [len(x) + 1 for x in seq_list]
dp = np.zeros(shape=shape, dtype=np.int64)

bt = np.zeros(shape=shape, dtype=np.int64)

bt_ = {
    0: (0, 0, 0, -1),
    1: (0, 0, -1, 0),
    2: (0, 0, -1, -1),
    3: (0, -1, 0, 0),
    4: (0, -1, 0, -1),
    5: (0, -1, -1, 0),
    6: (0, -1, -1, -1),
    7: (-1, 0, 0, 0),
    8: (-1, 0, 0, -1),
    9: (-1, 0, -1, 0),
    10: (-1, 0, -1, -1),
    11: (-1, -1, 0, 0),
    12: (-1, -1, 0, -1),
    13: (-1, -1, -1, 0),
    14: (-1, -1, -1, -1), # match/mismatch
}

def s(sl, up, ind, t=False):
    l = ''
    for i in range(4):
        l += '-' if up[i] == 0 else sl[i][ind[i] + up[i]]

    score = 0
    for i in range(4):
        for j in range(i + 1, 4):
            if l[i] != l[j]:
                score -= 1
    if t:
        return l, score
    return score

ln = [len(x) + 1 for x in seq_list]
index = [1, 1, 1, 1]

for index, x in np.ndenumerate(dp):
    if 0 in index:
        dp[index] = 0 if sum(index) == 0 else -100000
        if sum(index) == 0:
            dp[index] = 0
    else:
        bt[tuple(index)], dp[tuple(index)] = max([(k, dp[tuple([index[i] + v[i] for i in range(4)])] + s(seq_list, v, index)) for k, v in bt_.items()], key=lambda xx: xx[1])
        
alignments = ['' for i in seq_list]
index = [len(x) for x in seq_list]

print(dp[tuple(index)])

while any(index):
    update = bt_[bt[tuple(index)]]
    for x in range(4):
        if update[x] == -1:
            alignments[x] += seq_list[x][index[x] - 1]
        else:
            alignments[x] += '-'
        index[x] += update[x]


alignments = [x[::-1].replace('$', '') for x in alignments]

print('\n'.join(alignments))