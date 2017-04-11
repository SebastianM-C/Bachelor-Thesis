#!/usr/bin/env python

import numpy as np
from scipy import special as spec
from scipy import linalg
from timeit import default_timer as timer

from tools import cd


def elem(m, n, k, l):
    """
    Compute the matrix element on harmonic oscillator states
             + k    l
     < m | (a )  (a)  | n >
    """
    if m != n - l + k:
        return 0
    if n < l:
        return 0
    return 1 / spec.factorial(n - l) * \
        np.sqrt(spec.factorial(n) * spec.factorial(n - l + k))


a, b, d, n = 1, 0.2, 0.4, 60
cd(b, d, n)
nn = int(n * (n + 1) / 2)
H = np.empty((nn, nn))

index = [(n1, n2) for n1 in range(n) for n2 in range(n - n1)]

start = timer()
# Compute the Hamiltonian matrix elements
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        m1 = index[i][0]
        m2 = index[i][1]
        n1 = index[j][0]
        n2 = index[j][1]
        H[i][j] = a * (elem(m1, n1, 1, 1) * elem(m2, n2, 0, 0) +
                       elem(m1, n1, 0, 0) * elem(m2, n2, 1, 1)) \
            + 0.25 * b * (3 * elem(m1, n1, 1, 0) * elem(m2, n2, 2, 0)
                          + 3 * elem(m1, n1, 0, 1) * elem(m2, n2, 0, 2)
                          - elem(m1, n1, 3, 0) * elem(m2, n2, 0, 0)
                          - elem(m1, n1, 0, 3) * elem(m2, n2, 0, 0)) \
            + 0.75 * b * (elem(m1, n1, 0, 1) * elem(m2, n2, 2, 0)
                          + elem(m1, n1, 1, 0) * elem(m2, n2, 0, 2)
                          - elem(m1, n1, 1, 2) * elem(m2, n2, 0, 0)
                          - elem(m1, n1, 2, 1) * elem(m2, n2, 0, 0)
                          + 2 * elem(m1, n1, 0, 1) * elem(m2, n2, 1, 1)
                          + 2 * elem(m1, n1, 1, 0) * elem(m2, n2, 1, 1)) \
            + 0.375 * d * (elem(m1, n1, 2, 2) * elem(m2, n2, 0, 0)
                           + elem(m1, n1, 0, 0) * elem(m2, n2, 2, 2)) \
            + 0.125 * d * (elem(m1, n1, 2, 0) * elem(m2, n2, 0, 2)
                           + elem(m1, n1, 0, 2) * elem(m2, n2, 2, 0)) \
            + 0.500 * d * elem(m1, n1, 1, 1) * elem(m2, n2, 1, 1) \
            + 0.250 * d * (elem(m1, n1, 1, 3) * elem(m2, n2, 0, 0)
                           + elem(m1, n1, 3, 1) * elem(m2, n2, 0, 0)
                           + elem(m1, n1, 0, 0) * elem(m2, n2, 1, 3)
                           + elem(m1, n1, 0, 0) * elem(m2, n2, 3, 1)
                           + elem(m1, n1, 0, 2) * elem(m2, n2, 1, 1)
                           + elem(m1, n1, 2, 0) * elem(m2, n2, 1, 1)
                           + elem(m1, n1, 1, 1) * elem(m2, n2, 0, 2)
                           + elem(m1, n1, 1, 1) * elem(m2, n2, 2, 0)) \
            + 0.0625 * d * (elem(m1, n1, 4, 0) * elem(m2, n2, 0, 0)
                            + elem(m1, n1, 0, 4) * elem(m2, n2, 0, 0)
                            + elem(m1, n1, 0, 0) * elem(m2, n2, 4, 0)
                            + elem(m1, n1, 0, 0) * elem(m2, n2, 0, 4)
                            + 2 * elem(m1, n1, 2, 0) * elem(m2, n2, 2, 0)
                            + 2 * elem(m1, n1, 0, 2) * elem(m2, n2, 0, 2))

end = timer()
print(end - start)
np.savez_compressed('hamilt.npz', H=H)
# print(H)
# linalg.eigh(H)
