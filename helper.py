from typing import Generator
import numpy as np


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

def lps(pattern):
    i, j = 1, 0
    lps = [0 for x in range(len(pattern))]
    lps[0] = 0 # I know that's already 0 :))

    while i < len(pattern):
        if pattern[i] == pattern[j]:
            j += 1
            lps[i] = j
            i += 1
        else:
            if j != 0:
                j = lps[j - 1]
            else:
                lps[i] = 0
                i += 1

    return lps

def kmp(ref, sub) -> Generator:
    lps_ = lps(sub)
    i, j = 0, 0
    while i < len(ref):
        if ref[i] == sub[j]:
            i += 1
            j += 1
            if j == len(sub):
                yield i - j + 1
                j = lps_[j - 1]
        else:
            if j != 0:
                j = lps_[j - 1]
            else:
                i += 1

def lcs(x, y):
    dp = np.ndarray(shape=(len(x) + 1, len(y) + 1), dtype=np.str0)
    x_, y_ = len(x), len(y)
    for i in range(x_):
        dp[i, 0] = ''
    for i in range(y_):
        dp[0, i] = ''
    for i in range(1, x_ + 1):
        for j in range(1, y_ + 1):
            if x[i - 1] == y[j - 1]:
                dp[i, j] = dp[i - 1, j - 1] + x[i - 1]
            else:
                dp[i, j] = dp[i - 1, j] if len(dp[i - 1, j]) >= len(dp[i, j - 1]) else dp[i, j - 1]
    print(dp)
    return dp[x_, y_]
            