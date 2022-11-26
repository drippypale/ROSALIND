import os
import pathlib

with open(os.path.join(pathlib.Path(__file__).parent.resolve(), 'rna_codon.txt'), 'r') as f:
    lines = f.read().splitlines()
    RNA_CODON = {c.split(' ')[0]: c.split(' ')[1] for c in lines}