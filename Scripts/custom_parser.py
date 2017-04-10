import argparse
import numpy as np


def parse():
    """Parse arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', type=np.float64, nargs='+', default=[0.2],
                        help='Hamiltonian B parameter')
    parser.add_argument('-d', type=np.float64, nargs='+', default=[0.4],
                        help='Hamiltonian D parameter')
    parser.add_argument('-n', type=np.int64, nargs='+', default=[4],
                        help='diagonalisation basis size')
    parser.add_argument('-sc', action='store_true', default=False,
                        help='Specify whether to use SciPy or not for \
                      diagonalisation')
    # args = parser.parse_args(input().split())
    args = parser.parse_args()

    # Hamiltonian parameters
    B = args.b
    D = args.d
    N = args.n
    use_sc = args.sc    # Optionally use SciPy for eigenvalues and eigenvectors

    return B, D, N, use_sc
