'''
Drippypale
drippypale@gmail.com

Problem:
    name: Complementing a Strand of DNA
    link: https://rosalind.info/problems/revc/

'''
dna_s: str

with open('revc.in', 'r') as f:
    dna_s = f.read().splitlines()[0]

result = ''

complement = {
    'A':'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

for n in dna_s[::-1]:
    result += complement[n]
print(result)