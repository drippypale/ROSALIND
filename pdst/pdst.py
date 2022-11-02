'''
Drippypale
drippypale@gmail.com

Problem:
    name: Creating a Distance Matrix
    link: https://rosalind.info/problems/pdst/

'''

from ast import Num
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

sequence_list = list(FASTA.read_seq_list('pdst.in').values())
n = len(sequence_list) # number of sequences
l = len(sequence_list[0]) # as the problem mentioned, the sequences are all the same length so we store the length in l

import numpy as np

pdistance_matrix = np.zeros((n, n), dtype=np.float64)

for i in range(n):
    for j in range(i):
        t = 0 # stores the number of different nucleotides
        for c in range(l):
            if sequence_list[i][c] != sequence_list[j][c]:
                t += 1
        pdistance_matrix[i][j] = t / l
        pdistance_matrix[j][i] = pdistance_matrix[i][j] # the matrix is a diagonal matrix

for i in range(n):
    print(' '.join([f'{x:.7f}' for x in pdistance_matrix[i]]))
