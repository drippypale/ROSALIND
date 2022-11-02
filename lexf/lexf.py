'''
Drippypale
drippypale@gmail.com

Problem:
    name: Enumerating k-mers Lexicographically
    link: https://rosalind.info/problems/lexf/

'''

from collections import deque

alphabet: list
n: int

with open('lexf.in', 'r') as f:
    alphabet = f.readline()[:-1].split(' ')
    alphabet.sort()
    n = int(f.readline()[:-1])


strings = deque(alphabet) #linked list O(1) on delete/insert

for i in range(1, n):
    for j in range(len(strings)):
        s = strings.popleft()
        for c in alphabet:
            strings.append(s + c)

print('\n'.join(strings))