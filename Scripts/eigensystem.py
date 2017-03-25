#!/usr/bin/env python
import os
import numpy as np
from scipy import linalg
from timeit import default_timer as timer
import matplotlib.pyplot as plt


def get(use_sc, return_H=False):
    """Return the eigenvalues and the eigenvectors (along with the state index
    and max coefficient index) and optionally the Hamiltonian"""
    # Load files
    H = np.loadtxt("hamilt.out")    # transposed Hamiltonian
    index = np.loadtxt("index.out", dtype=int)
    print("Now in: ", os.getcwd())
    H = np.transpose(H)
    # Get eigenvalues and eigenvectors
    if use_sc:
        E, eigenvectors = linalg.eigh(H)
    else:
        E = np.loadtxt("hamilt.dat", usecols=1, unpack=True)    # energy levels
        eigenvectors = np.loadtxt("eigenvectors.out", unpack=True)

    eigenvectors = np.transpose(eigenvectors)  # each eigenvector is on one row

    # max coefficient in eigenvector
    c_max = np.empty_like(eigenvectors[0], dtype=int)

    # The index of the largest coefficient
    for i in range(eigenvectors[0].size):
        c_max[i] = np.argmax(np.abs(eigenvectors[i]))

    if return_H:
        return E, eigenvectors, index, c_max, H
    return E, index[c_max]


def eigv(use_sc):
    """Return the eigenvalues"""
    H = np.loadtxt("hamilt.out")    # transposed Hamiltonian
    H = np.transpose(H)
    # Get eigenvalues and eigenvectors
    if use_sc:
        E = linalg.eigvalsh(H)
    else:
        E = np.loadtxt("hamilt.dat", usecols=1, unpack=True)    # energy levels
    return E


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
    epsilon = 1e-3 if not use_sc else 1e-5
    levels = np.split(E, np.where(np.diff(E) > epsilon)[0] + 1)

    plt.hist(np.diff(E),
             bins=[0, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3,
                   5e-3, 1e-2, 5e-2, 0.1, 0.5, 1]
             )
    plt.xscale('log', nonposy='clip')
    plt.savefig('np.diff(E)' + ('_sc.png' if use_sc else '.png'))
    plt.close()

    k = 0
    for i in range(len(levels)):
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
