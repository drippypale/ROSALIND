'''
Drippypale
drippypale@gmail.com

Problem:
    name: Enumerating Gene Orders
    link: https://rosalind.info/problems/perm/

'''

import math
import itertools

n: int
with open('perm.in', 'r') as f:
    n = int(f.read().splitlines()[0])

print(math.factorial(n))

for i in itertools.permutations(range(1, n + 1)):
    print(' '.join([str(x) for x in i]))
