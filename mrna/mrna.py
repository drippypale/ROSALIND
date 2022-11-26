'''
Drippypale
drippypale@gmail.com

Problem:
    name: Inferring mRNA from Protein
    link: https://rosalind.info/problems/mrna/

'''

from utils.rna_codon import RNA_CODON # a dict of rna codons according to the problem page

codon_rev = dict()
for v in RNA_CODON.values():
    if codon_rev.get(v, False):
        codon_rev[v] += 1
    else:
        codon_rev[v] = 1

with open('mrna.in', 'r') as f:
    protein = f.read().splitlines()[0]

result = 1
for p in protein:
    result *= codon_rev[p]
    if result >= 1_000_000:
        result -= 1_000_000

result *= codon_rev['Stop']
print(result % 1_000_000)
