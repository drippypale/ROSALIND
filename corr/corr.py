'''
Drippypale
drippypale@gmail.com

Problem:
    name: Error Correction in Reads
    link: https://rosalind.info/problems/corr/

'''

from utils.helper import FASTA
complement = {
    'A':'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

from collections import deque
# seq_list = [[seq, False] for seq in FASTA.read_seq_list('corr.in').values()]
seq_list = list(FASTA.read_seq_list('corr.in').values())
n = len(seq_list)

revc = lambda seq: ''.join([complement[c] for c in seq[::-1]])
distance = lambda seq1, seq2: sum([0 if seq1[i] == seq2[i] else 1 for i in range(len(seq1))])

seq_list_revc = [revc(seq) for seq in seq_list]

def corr(seq, seq_list, seq_list_revc):
    for seq_ in seq_list:
        if distance(seq_, seq) == 1:
            return seq_
    for seq_ in seq_list_revc:
        if distance(seq_, seq) == 1:
            return seq_

solid = [seq_ for seq_ in seq_list if seq_list.count(seq_) > 1 or seq_ in seq_list_revc]
solid_revc = [revc(seq_) for seq_ in solid]

corrections = [f'{seq}->{corr(seq, solid, solid_revc)}' for seq in seq_list if seq not in solid]

print('\n'.join(corrections))


