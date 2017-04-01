#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import diff
import eigensystem
from tools import cd
from tools import clean_dir
from custom_parser import parse


def main(b, d, n, use_sc):
    """Create bar plots for eigenvectors"""
    cd(b, d, n)

    E, eigenvectors, ket = eigensystem.eigv(use_sc)
    stable_levels = np.load('cache.npy')    # get cached stable levels
    # Select only stable levels
    eigenvectors = eigenvectors[:stable_levels]
    ket = ket[:stable_levels]
    # Get the index array that sorts the eigenvectors with respect to n1 and n2
    sort_idx = np.argsort(ket.view('<i4, <i4'),
                          order=['f0', 'f1'], axis=0).reshape(ket.shape[0],)
    # Sort the eigenvectors and the eigenvalues
    ket = ket[sort_idx]
    eigenvectors = eigenvectors[sort_idx]
    E = E[sort_idx]

    clean_dir('eigenvectors')
    for i in range(eigenvectors.shape[0]):
        plt.bar(range(eigenvectors.shape[1]), np.abs(eigenvectors[i]),
                label='E = ' + str(E[i]))
        plt.legend()
        plt.savefig('eigenvectors/ket' + str(ket[i][0]) + ' ' + str(ket[i][1]),
                    dpi=900)
        # plt.show()
        plt.close()


if __name__ == '__main__':
    B, D, N, use_sc = parse()
    for b in B:
        for d in D:
            for n in N:
                main(b, d, n, use_sc)
