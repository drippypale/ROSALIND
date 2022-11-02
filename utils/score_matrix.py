from typing import Tuple
import numpy as np

import os
import pathlib

def _read_score_matrix(path) -> Tuple[dict, np.ndarray]:
    alphabet_index: dict
    matrix: np.ndarray

    with open(path, 'r') as f:
        lines = f.read().splitlines()
        alphabet_index = {c: i for i, c in enumerate(lines[0].split(' '))}
        n = len(alphabet_index)
        matrix = np.zeros(shape=(n, n), dtype=np.int64)
        for i in range(1, n + 1):
            for j, j_ in enumerate(lines[i].split(' ')):
                matrix[i - 1, j] = int(j_)
    return alphabet_index, matrix


BLOSUM62            = _read_score_matrix(os.path.join(pathlib.Path(__file__).parent.resolve(), 'score_matrices/BLOSUM62'))
PAM250              = _read_score_matrix(os.path.join(pathlib.Path(__file__).parent.resolve(), 'score_matrices/PAM250'))
