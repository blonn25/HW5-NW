import numpy as np


seqA = "AGCT"
seqB = "AGTT"
len_A = len(seqA)
len_B = len(seqB)
gap_open = -2
gap_extend = -1

_align_matrix = np.zeros((len_A + 1, len_B + 1))
_align_matrix[0, 0] = 0
for i in range(1, len_A + 1):
    _align_matrix[i, 0] = gap_open + (i-1) * gap_extend
for j in range(1, len_B + 1):
    _align_matrix[0, j] = gap_open + (j-1) * gap_extend

print(_align_matrix[(1,0)])