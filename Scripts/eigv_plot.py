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
    # Get states
    state, ket, index = eigensystem.get_state(use_sc)
    # Get the index array that sorts the eigenvector coefficients
    # such that n1 + n2 is incerasing
    sort_idx = np.argsort(index.sum(axis=1))
    # Sort the eigenvector coefficients
    state['eigvec'] = state['eigvec'][:, sort_idx]
    # stable_levels = np.load('cache.npy')    # get cached stable levels
    stable_levels = 10
    # Select the states corresponding to stable levels
    state = state[:stable_levels]
    ket = ket[:stable_levels]

    # Get irreductible representation index
    ir_reps = eigensystem.levels(state['E'], ket, use_sc)
    # Build irreductible representation string
    reps = 'reuns', 'reuna', 'rebde'
    ir_str = [reps[i] for i in ir_reps]

    # Compute energy plot parameters
    k = index[sort_idx].sum(axis=1)         # n1 + n2
    w = 1 / (k + 1) - 0.01                  # bar widths
    r = [j for i in range(n + 1) for j in range(i)]
    x = k - 0.5 + r / (k + 1)               # positions
    x = x[:stable_levels]
    w = w[:stable_levels]

    clean_dir('eigenvectors')
    minor_ticks = np.arange(0, 10, 0.5)
    for i in range(state.shape[0]):
        # Index plot
        plt.bar(range(10), np.abs(state['eigvec'][i, :10]), label='E = ' +
                str(state['E'][i]) + '\n' + '$\\left|' +
                str(ket[i][0]) + str(ket[i][1]) + '\\right\\rangle$\t' +
                ir_str[i])
        plt.legend()
        plt.xticks(range(10),
                   ['$c_{' + str(index[sort_idx][i][0]) +
                    str(index[sort_idx][i][1]) + '}$'
                    for i in range(index.shape[0])],
                   )
        plt.ylabel('$C_{ij}$')
        plt.savefig('eigenvectors/ket_i' + str(ket[i][0]) + ' '
                    + str(ket[i][1]), dpi=300)
        # plt.show()
        plt.close()
        # Energy plot
        plt.bar(x, np.abs(state['eigvec'][i, :10]), width=w, align='edge',
                label='E = ' + str(state['E'][i]) + '\n' + '$\\left|' +
                str(ket[i][0]) + str(ket[i][1]) + '\\right\\rangle$\t' +
                ir_str[i])
        plt.legend()
        plt.xticks(range(10), ['$E_{' + str(i) + '}$' for i in range(10)])
        plt.axes().set_xticks(minor_ticks, minor=True)
        plt.ylabel('$C_{ij}$')
        plt.savefig('eigenvectors/ket_e' + str(ket[i][0]) + ' ' +
                    str(ket[i][1]), dpi=300)
        # plt.show()
        plt.close()


if __name__ == '__main__':
    B, D, N, use_sc = parse()
    for b in B:
        for d in D:
            for n in N:
                main(b, d, n, use_sc)
