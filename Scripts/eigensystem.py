#!/usr/bin/env python

import numpy as np
from scipy import linalg
from scipy.io import FortranFile
from timeit import default_timer as timer
import matplotlib.pyplot as plt

from tools import get_input


def readH(format):
    """Read the Hamiltonian using the given format"""
    if format == 'fortran_bin':
        _, _, n = get_input()
        nn = int(n * (n + 1) / 2)
        hamilt = FortranFile('hamilt.bin', 'r')
        H = np.empty((nn, nn))
        for i in range(nn):
            H[i] = hamilt.read_reals(dtype='float32').reshape(nn)
        return H
    if format == 'text':
        H = np.loadtxt("hamilt.out")
        return H


def get(use_sc, return_H=False):
    """Return the eigenvalues and the eigenvectors (along with the state index
    and max coefficient index) and optionally the Hamiltonian"""
    # Load files
    H = readH('fortran_bin').T    # read transposed Hamiltonian
    index = np.loadtxt("index.out", dtype=int)
    # Get eigenvalues and eigenvectors
    if use_sc:
        E, eigenvectors = linalg.eigh(H)
    else:
        E = np.loadtxt("hamilt.dat", usecols=1, unpack=True)    # energy levels
        eigenvectors = np.loadtxt("eigenvectors.out", unpack=True)

    eigenvectors = np.transpose(eigenvectors)  # each eigenvector is on one row

    # max coefficient in eigenvector
    c_max = np.empty(eigenvectors.shape[0], dtype=int)

    # The index of the largest coefficient
    for i in range(eigenvectors.shape[0]):
        c_max[i] = np.argmax(np.abs(eigenvectors[i]))

    if return_H:
        return E, eigenvectors, index, c_max, H
    return E, index[c_max]


def get_state(use_sc):
    """Return the eigenvalues and the eigenvectors"""
    # Load files
    H = readH('fortran_bin').T    # read transposed Hamiltonian
    index = np.loadtxt("index.out", dtype=int)
    # Define the state array
    state = np.zeros(H.shape[0],
                     [('E', np.float64), ('eigvec', np.float64, H.shape[0])])
    # Get eigenvalues and eigenvectors
    if use_sc:
        state['E'], state['eigvec'] = linalg.eigh(H)
    else:
        state['E'] = np.loadtxt("hamilt.dat", usecols=1, unpack=True)
        state['eigvec'] = np.loadtxt("eigenvectors.out", unpack=True)

    # each eigenvector is on one row
    state['eigvec'] = np.transpose(state['eigvec'])

    # max coefficient in eigenvector
    c_max = np.empty_like(state['eigvec'][0], dtype=int)

    # The index of the largest coefficient
    for i in range(state['eigvec'][0].size):
        c_max[i] = np.argmax(np.abs(state['eigvec'][i]))

    return state, index[c_max], index


def levels(E, ket, use_sc, colors=''):
    """Return the degenerate subspace index and optionally the colormap"""
    # Irreductible representations
    # 0 - unidimensional symmetric representation (reuns)
    # 1 - unidimensional anti-symmetric representation (reuna)
    # 2 - bidimensional representation (rebde)
    ir_reps = np.zeros([E.size], dtype=np.uint8)
    return_colors = len(colors)
    if return_colors:
        colormap = [''] * E.size   # colors used

    # Group energy levels such that a level contains all the eigenvalues with
    # the same value
    # epsilon = 1e2   # if use_sc else 5e-2
    epsilon = 1e-5
    delta = np.diff(E)
    # delta_left = np.pad(delta, (1, 0), 'constant', constant_values=1)[:-1]
    # delta_right = np.pad(delta, (0, 1), 'constant', constant_values=1)[1:]
    # delta_avg = (delta_left + delta_right) / 2
    # delta_avg = np.sqrt(delta_left * delta_right)

    # levels = np.split(E, np.where(delta_avg / delta < epsilon)[0] + 1)
    levels = np.split(E, np.where(delta > epsilon)[0] + 1)

    # plt.hist(delta_avg / delta,
    #          bins=np.pad(np.geomspace(1e-9, 1e9, 19), (1, 0), mode='constant'),
    #          label='$\\frac{\\sqrt{\\Delta_l\\,\\Delta_r}}{\\Delta}$'
    #          )
    # plt.legend()
    # plt.xscale('log', nonposy='clip')
    # plt.savefig('delta_r' + ('_sc.png' if use_sc else '.png'))
    # plt.show()
    # plt.close()
    # Energy difference (between two consecutive levels) histogram
    plt.hist(delta,
             bins=np.pad(np.geomspace(1e-9, 1e3, 13), (1, 0), mode='constant'),
             label='$\\Delta = E_{n+1} - E_n$'
             )
    plt.legend()
    plt.xscale('log')
    plt.savefig('hist_delta' + ('_sc.png' if use_sc else '.png'))
    # plt.show()
    plt.close()
    # Energy difference bar plot
    plt.figure(figsize=(20, 4))
    plt.bar(range(1000), delta[:1000], label='$\\Delta = E_{n+1} - E_n$')
    plt.axhline(y=epsilon)
    plt.legend()
    plt.yscale('log', nonposy='clip')
    # plt.xscale('log')
    plt.savefig('bar_delta' + ('_sc.png' if use_sc else '.png'), dpi=600,
                bbox_inches='tight')
    plt.show()
    plt.close()

    k = 0
    for i in range(len(levels)):
        # Check for bidimensional representation selection problems
        if levels[i].size > 2:
            print('Warning: bidimensional representation selection',
                  'problem: size: ', levels[i].size, 'at', i,
                  '\nenergy: ', levels[i], '\ndelta: ', np.diff(levels[i]))
        for j in range(levels[i].size):
            if return_colors:
                colormap[i + j + k] = colors[i % len(colors)]
            if levels[i].size > 1:  # degenerate subspace -> rebde
                ir_reps[i + j + k] = 2
            else:
                if ket[i][1] % 2:   # n2 odd -> reuna
                    ir_reps[i + j + k] = 1
                else:               # n2 even -> reuns
                    ir_reps[i + j + k] = 0
        k += levels[i].size - 1

    if return_colors:
        return ir_reps, colormap
    return ir_reps
