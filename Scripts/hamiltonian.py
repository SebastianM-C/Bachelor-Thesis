#!/usr/bin/env python

import numpy as np
# import sympy as sp
from scipy import special
from timeit import default_timer as timer

from tools import cd
from custom_parser import parse


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
    # Use sympy to simplify the expression
    # F = sp.factorial
    # h = 1 / F(n - l) * sp.sqrt(F(n) * F(n - l + k))
    # return h.simplify()
    F = special.factorial

    return 1 / np.float64(F(n - l, exact=True)) * \
        np.sqrt(np.float64(F(n, exact=True)) *
                np.float64(F(n - l + k, exact=True)))


def main(a, b, d, n):
    n = int(n)
    cd(b, d, n)
    nn = int(n * (n + 1) / 2)
    H = np.empty((nn, nn))

    index = [(n1, n2) for n1 in range(n) for n2 in range(n - n1)]

    np.seterr(over='raise')
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
    print('hamilt: ', end - start)
    np.savez_compressed('hamilt.npz', H=H)
    return H
    # print(H)


if __name__ == '__main__':
    B, D, N, _ = parse()
    for b in B:
        for d in D:
            for n in N:
                main(1, b, d, n)
