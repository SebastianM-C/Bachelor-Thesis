# !/usr/bin/env python

import os
import numpy as np
import eigensystem

import tools
from plots import bar_plot, histogram


def relSpacing(E):
    deltaE = np.diff(E)
    avgSpacing = (E[-1] - E[0]) / E.size
    return deltaE / avgSpacing


def difference(E1, b, d, n, delta_n, use_sc=False, ir_reps1=np.empty(0)):
    """Return the differences between the energy levels and of the
    structure of the irreductible representations (optional) for the given
    diagonalization basis and n + delta_n"""
    # Load the energy levels in the second basis
    os.chdir("../B" + str(b) + " D" + str(d) + " N" + str(n + delta_n))
    E2, ket2 = eigensystem.get(use_sc, return_ket=True)
    if ir_reps1.size:
        ir_reps2 = eigensystem.levels(E2, ket2, use_sc)

    if E2.size > E1.size:
        E_diff = (E2[:E1.size] - E1) / E1
        if ir_reps1.size:
            ir_diff = ir_reps2[:ir_reps1.size] - ir_reps1
    else:
        E_diff = (E1[:E2.size] - E2) / E2
        if ir_reps1.size:
            ir_diff = ir_reps1[:ir_reps2.size] - ir_reps2

    os.chdir("../B" + str(b) + " D" + str(d) + " N" + str(n))
    if ir_reps1.size == 0:
        return np.abs(E_diff)
    return np.abs(E_diff), ir_diff


def stable(E1, b, d, n, use_sc, delta_n, epsilon, ir_reps=np.empty(0)):
    """Return the number of stable levels"""
    if ir_reps.size == 0:
        E_diff = difference(E1, b, d, n, delta_n, use_sc)
    else:
        E_diff, ir_diff = difference(E1, b, d, n, delta_n, use_sc, ir_reps)
    # Energy difference (between two diagonalization bases) histogram
    histogram(E_diff, label='B' + str(b) + ' D' + str(d) + ' N' + str(n),
              bins=np.pad(np.geomspace(1e-14, 1e-2, 13), (1, 0),
                          mode='constant'), xscale='log',
              fname='hist_E_diff' + ('_sc.png' if use_sc else '.png')
              )
    # Energy difference bar plot
    bar_plot(E_diff[E_diff < 0.01],
             label='B' + str(b) + ' D' + str(d) + ' N' + str(n),
             figsize=(20, 4), axhline_y=epsilon, yscale='log', dpi=600,
             fname='bar_E_diff' + ('_sc.png' if use_sc else '.png'),
             bbox_inches='tight')

    print("E_diff > epsilon :",
          np.where(E_diff > epsilon)[0][:5], "\nstability epsilon: ", epsilon)
    if ir_reps.size:
        print("ir_diff: ", np.where(ir_diff > 0)[0][:5])
    # Cache the result
    np.where(E_diff > epsilon)[0][1].tofile('stable.txt', sep=' ')
    return np.where(E_diff > epsilon)[0][1]
    # return np.where(ir_diff > 0)[0][1]
