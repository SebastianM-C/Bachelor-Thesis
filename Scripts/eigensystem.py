#!/usr/bin/env python

import numpy as np
from scipy import linalg
from timeit import default_timer as timer
from os.path import isfile

from tools import get_input
from hamiltonian import main as hamiltonian
from plots import bar_plot, histogram


def readH(format):
    """Read the Hamiltonian using the given format"""
    if format == 'npz':
        if isfile('hamilt.npz'):
            hamilt = np.load('hamilt.npz')
            return hamilt['H']
        else:
            # Fallback to Fortran binary
            format = 'fortran_bin'
            # print('Hamiltonian file not found. Computing again.')
            # b, d, n = get_input()
            # return hamiltonian(1, b, d, n)
    if format == 'fortran_bin':
        _, _, n = get_input()
        nn = int(n * (n + 1) / 2)
        H = np.empty((nn, nn))
        with open('hamilt.bin', 'rb') as h:
            for i in range(nn):
                H[i] = np.fromfile(h, dtype='float64', count=nn).reshape(nn)
        return H
    if format == 'text':
        H = np.loadtxt("hamilt.out")
        return H.T


def get(use_sc, return_eigv=False, return_ket=False, return_index=False,
        return_cmax=False, return_H=False):
    """Return the eigenvalues and optionally the eigenvectors,
    the number operator form of the states(ket), the state index of the states,
    the max coefficient index and the Hamiltonian"""
    # Load files
    H = readH('npz')    # read the Hamiltonian
    # Save to npz to save sapce
    if not isfile('hamilt.npz'):
        np.savez_compressed('hamilt.npz', H=H)
    b, d, n = get_input()
    n = int(n)
    index = np.array([(n1, n2) for n1 in range(n) for n2 in range(n - n1)])
    # Get eigenvalues and eigenvectors
    if isfile('eigensystem.npz'):
        print('Used cached result for: B =', b, ' D =', d, ' N =', n)
        eigensystem = np.load('eigensystem.npz')
        E = eigensystem['E']
        eigenvectors = eigensystem['eigenvectors']
    else:
        start = timer()
        E, eigenvectors = linalg.eigh(H)
        end = timer()
        print('Diagonalisation for N =', n, ':', end - start, 'seconds')
        # Save the results
        np.savez_compressed('eigensystem.npz', E=E, eigenvectors=eigenvectors)
    # else:
    #     E = np.loadtxt("hamilt.dat", usecols=1, unpack=True)    # energy levels
    #     eigenvectors = np.loadtxt("eigenvectors.out", unpack=True)

    eigenvectors = np.transpose(eigenvectors)  # each eigenvector is on one row

    # max coefficient in eigenvector
    c_max = np.empty(eigenvectors.shape[0], dtype=int)

    # The index of the largest coefficient
    for i in range(eigenvectors.shape[0]):
        c_max[i] = np.argmax(np.abs(eigenvectors[i]))

    results = (E, )
    if return_eigv:
        results = results + (eigenvectors, )
    if return_ket:
        results = results + (index[c_max], )
    if return_index:
        results = results + (index, )
    if return_cmax:
        results = results + (c_max, )
    if return_H:
        results = results + (H, )
    return results


def levels(E, ket, use_sc, epsilon=1e-8, colors=''):
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
    delta = np.diff(E)
    avgSpacing = (E[-1] - E[0]) / E.size
    relsp = delta / avgSpacing
    print('levels epsilon: ', epsilon)
    print('avgSpacing: ', avgSpacing)

    levels = np.split(E, np.where(relsp > epsilon)[0] + 1)

    # Energy difference (between two consecutive levels) histogram
    histogram(delta, label='$\\Delta E$', xscale='log',
              bins=np.pad(np.logspace(-15, 1, 17), (1, 0),
                          mode='constant'), fname='hist_delta.png')
    # Relative spacing histogram
    histogram(relsp, label='$N \\frac{\\Delta E}{E_n - E_0}$', xscale='log',
              bins=np.pad(np.logspace(-13, 1, 15), (1, 0),
                          mode='constant'), fname='hist_relsp.png', xlabel='S')
    # Energy difference bar plot
    bar_plot(delta, figsize=(20, 4), label='$\\Delta E$', yscale='log',
             fname='bar_delta' + ('_sc.png' if use_sc else '.png'), dpi=600,
             bbox_inches='tight')
    # Relative spacing bar plot
    bar_plot(relsp, figsize=(20, 4), label='$N \\frac{\\Delta E}{E_n - E_0}$',
             yscale='log', fname='relsp' + ('_sc.png' if use_sc else '.png'),
             dpi=600, axhline_y=epsilon, bbox_inches='tight', ylabel='S',
             show=False)

    k = 0
    for i in range(len(levels)):
        # Check for bidimensional representation selection problems
        if levels[i].size > 2:
            print('Warning: bidimensional representation selection',
                  'problem: size: ', levels[i].size, 'at', i,
                  '\nenergy: ', levels[i], '\ndelta: ', np.diff(levels[i]),
                  '\nrelsp: ', np.diff(levels[i]) / avgSpacing)
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
