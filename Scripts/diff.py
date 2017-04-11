#!/usr/bin/env python

import os
import numpy as np
import eigensystem
import matplotlib.pyplot as plt

import tools


def difference(E1, ir_reps1=None, b, d, n, delta_n, use_sc=False):
    """Return the differences between the energy levels and of the
    structure of the irreductible representations (optional) for the given
    diagonalization basis and n + delta_n"""
    # Load the energy levels in the second basis
    os.chdir("../B" + str(b) + " D" + str(d) + " N" + str(n + delta_n))
    E2, ket2 = eigensystem.get(use_sc)
    if ir_reps is not None:
        ir_reps2 = eigensystem.levels(E2, ket2, use_sc)

    if E2.size > E1.size:
        E_diff = (E2[:E1.size] - E1) / E1
        if ir_reps is not None:
            ir_diff = ir_reps2[:ir_reps1.size] - ir_reps1
    else:
        E_diff = (E1[:E2.size] - E2) / E2
        if ir_reps is not None:
            ir_diff = ir_reps1[:ir_reps2.size] - ir_reps2

    os.chdir("../B" + str(b) + " D" + str(d) + " N" + str(n))
    if ir_reps is None:
        return np.abs(E_diff)
    return np.abs(E_diff), ir_diff


def stable(E1, ir_reps=None, b, d, n, use_sc, delta_n):
    """Return the number of stable levels"""
    epsilon = 1e-8 if use_sc else 5e-2
    if ir_reps is None:
        E_diff = difference(E1, b, d, n, delta_n, use_sc)
    else:
        E_diff, ir_diff = difference(E1, ir_reps, b, d, n, delta_n, use_sc)
    # Energy difference (between two diagonalization bases) histogram
    plt.hist(E_diff,
             bins=np.pad(np.geomspace(1e-14, 1e-2, 13), (1, 0),
                         mode='constant'),
             label='B' + str(b) + ' D' + str(d) + ' N' + str(n)
             )
    plt.xscale('log')
    plt.savefig('hist_E_diff' + ('_sc.png' if use_sc else '.png'))
    # plt.show()
    plt.close()
    # Energy difference bar plot
    plt.figure(figsize=(20, 4))
    plt.bar(range(1000), E_diff[:1000],
            label='B' + str(b) + ' D' + str(d) + ' N' + str(n))
    plt.axhline(y=epsilon)
    plt.legend()
    plt.yscale('log', nonposy='clip')
    plt.xscale('log')
    plt.savefig('bar_E_diff' + ('_sc.png' if use_sc else '.png'), dpi=600,
                bbox_inches='tight')
    # plt.show()
    plt.close()

    print("E_diff > epsilon :",
          np.where(E_diff > epsilon)[0][:5], "\nepsilon: ", epsilon)
    if ir_reps is not None:
        print("ir_diff: ", np.where(ir_diff > 0)[0][:5])
    # Cache the result
    np.where(E_diff > epsilon)[0][1].tofile('stable.txt', sep=' ')
    return np.where(E_diff > epsilon)[0][1]
    # return np.where(ir_diff > 0)[0][1]
