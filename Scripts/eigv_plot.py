#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from os.path import isfile

import eigensystem
from tools import cd
from tools import clean_dir
from custom_parser import parse


def no_signif_el(eigvec):
    """Determine the number of significant eigenvector elements based on the
    partial norm"""
    sq_sum = 0
    for i in range(eigvec.shape[0]):
        sq_sum += eigvec[i]**2
        if np.sqrt(sq_sum) >= 0.95:
            no_el = i + 3
            print('eigv_len: ', no_el, '\npartial norm:', np.sqrt(sq_sum))
            break
    if no_el < 10:
        no_el = 10
    return no_el


def index_plot(eigvec, eigv_len, label, index, sort_idx, fname, d_no):
    """Index plot"""
    plt.bar(range(eigv_len), np.abs(eigvec[:eigv_len]),
            label=label)
    plt.legend()
    plt.xticks(range(eigv_len),
               ['$c_{' + str(index[sort_idx][i][0]) +
                str(index[sort_idx][i][1]) + '}$'
                for i in range(index.shape[0])],
               )
    plt.ylabel('$C_{ij}$')
    # Prevent overwriting for repeated states
    if not isfile('eigenvectors/ket_i' + fname + '.png'):
        plt.savefig('eigenvectors/ket_i' + fname, dpi=400)
    else:
        d_no += 1
        plt.savefig('eigenvectors/ket_i' + fname + '_' + str(d_no),
                    dpi=400)
        print('|' + fname, '> is not unique')
    # plt.show()
    plt.close()


def energy_plot(eigvec, eigv_len, x, w, label, index, sort_idx, fname, d_no):
    """Energy plot"""
    plt.bar(x[:eigv_len], np.abs(eigvec[:eigv_len]),
            width=w[:eigv_len], align='edge', label=label)
    plt.legend()
    # Use the last position for the length of the axis
    x_len = int(x[:eigv_len][-1]) + 2
    minor_ticks = np.arange(0, x_len, 0.5)
    plt.xticks(range(x_len),
               ['$E_{' + str(j) + '}$' for j in range(x_len)])
    plt.axes().set_xticks(minor_ticks, minor=True)
    plt.ylabel('$C_{ij}$')
    # Prevent overwriting for repeated states
    if not isfile('eigenvectors/ket_e' + fname + '.png'):
        plt.savefig('eigenvectors/ket_e' + fname, dpi=400)
    else:
        plt.savefig('eigenvectors/ket_e' + fname + '_' + str(d_no),
                    dpi=400)
    # plt.show()
    plt.close()


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
    no_eigv = 40        # number of eigenvectors to plot
    # Select the states corresponding to stable levels
    state = state[:no_eigv]
    ket = ket[:no_eigv]
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

    d_no = 0   # number of duplicate states

    clean_dir('eigenvectors')
    for i in range(state.shape[0]):
        eigv_len = no_signif_el(state['eigvec'][i])
        minor_ticks = np.arange(0, eigv_len, 0.5)
        # Plot label
        label = 'E = ' + str(state['E'][i]) + '\n' + '$\\left|' + \
            str(ket[i][0]) + '\\,' + str(ket[i][1]) + '\\right\\rangle$\t' + \
            ir_str[i]
        fname = str(ket[i][0]) + ' ' + str(ket[i][1])   # filename

        index_plot(state['eigvec'][i], eigv_len, label, index, sort_idx, fname,
                   d_no)

        energy_plot(state['eigvec'][i], eigv_len, x, w, label, index, sort_idx,
                    fname, d_no)


if __name__ == '__main__':
    B, D, N, use_sc = parse()
    for b in B:
        for d in D:
            for n in N:
                main(b, d, n, use_sc)
