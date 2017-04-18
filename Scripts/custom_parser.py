import argparse
import numpy as np


def parse(advanced=False, select=False):
    """Parse arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', type=np.float64, nargs='+', default=[0.2],
                        help='Hamiltonian B parameter')
    parser.add_argument('-d', type=np.float64, nargs='+', default=[0.4],
                        help='Hamiltonian D parameter')
    parser.add_argument('-n', type=np.int64, nargs='+', default=[4],
                        help='Diagonalisation basis size')
    parser.add_argument('-sc', action='store_true', default=False,
                        help='Specify whether to use SciPy or not for' +
                        'diagonalisation')
    parser.add_argument('-dn', '--delta_n', type=np.int64, nargs=1,
                        default=20,
                        help='Size difference between diagonalization bases')
    parser.add_argument('-st_eps', '--stability_epsilon', type=np.float64,
                        default=1e-11,
                        help='Maximum difference between two stable levels' +
                        'when the diagonalization basis size is increased by' +
                        'delta_n')
    parser.add_argument('-l_eps', '--levels_epsilon', type=np.float64,
                        default=1e-8,
                        help='Minimum difference between two consecutive' +
                        'levels from one of the irreductible' +
                        'unidimensional representations')
    parser.add_argument('-r', '--reselect', action='store_false',
                        default=True, help='Specify whether to reselect the' +
                        'irreductible representations or not')

    # args = parser.parse_args(input().split())
    args = parser.parse_args()

    # Hamiltonian parameters
    B = args.b
    D = args.d
    N = args.n
    use_sc = args.sc    # Optionally use SciPy for eigenvalues and eigenvectors
    delta_n = args.delta_n
    st_epsilon = args.stability_epsilon
    lvl_epsilon = args.levels_epsilon
    reselect = args.reselect

    if advanced and select:
        return B, D, N, use_sc, delta_n, st_epsilon, lvl_epsilon, reselect
    if advanced:
        return B, D, N, use_sc, delta_n, st_epsilon, lvl_epsilon
    return B, D, N, use_sc
