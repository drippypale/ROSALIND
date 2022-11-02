'''
Drippypale
drippypale@gmail.com

Problem:
    name: Computing GC Content
    link: https://rosalind.info/problems/gc/

'''

lines = list()
with open('gc.in', 'r') as f:
    lines = f.read().splitlines()

max_ = (.0, '', '') # (GC_content, seq, tag)
t_seq: str = None
t_tag: str = None
for line in lines:
    if line[0] == '>':
        if t_seq is not None and (t_seq.count('C') + t_seq.count('G')) / len(t_seq) > max_[0]:
            max_ = ((t_seq.count('C') + t_seq.count('G')) / len(t_seq), t_seq, t_tag)
        t_seq, t_tag = '', line[1:]
    else:
        t_seq += line

if (t_seq.count('C') + t_seq.count('G')) / len(t_seq) > max_[0]:
        max_ = ((t_seq.count('C') + t_seq.count('G')) / len(t_seq), t_seq, t_tag)
    
print(f'{max_[2]}\n{max_[0]*100:.7f}')            
