#!/usr/bin/env python

import os
import numpy as np
import eigensystem
import matplotlib.pyplot as plt

import tools


def difference(E1, ir_reps1, b, d, n, delta_n, use_sc=False):
    """Return the differences between the energy levels and of the
    structure of the irreductible representations for the given
    diagonalization basis and n + delta_n"""
    # Load the energy levels in the second basis
    os.chdir("../B" + str(b) + " D" + str(d) + " N" + str(n + delta_n))
    E2, ket2 = eigensystem.get(use_sc)
    ir_reps2 = eigensystem.levels(E2, ket2, use_sc)

    if E2.size > E1.size:
        E_diff = (E2[:E1.size] - E1) / E1
        ir_diff = ir_reps2[:ir_reps1.size] - ir_reps1
    else:
        E_diff = (E1[:E2.size] - E2) / E2
        ir_diff = ir_reps1[:ir_reps2.size] - ir_reps2

    os.chdir("../B" + str(b) + " D" + str(d) + " N" + str(n))
    return np.abs(E_diff), ir_diff


def stable(E1, ir_reps, b, d, n, use_sc, delta_n):
    """Return the number of stable levels"""
    epsilon = 1e-2 if not use_sc else 1e-12
    # delta_n = 10 if not use_sc else 20
    E_diff, ir_diff = difference(E1, ir_reps, b, d, n, delta_n, use_sc)

    plt.hist(E_diff,
             bins=[0, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8,
                   1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
             label='B' + str(b) + ' D' + str(d) + ' N' + str(n)
             )
    plt.xscale('log', nonposy='clip')
    plt.savefig('E_diff' + ('_sc.png' if use_sc else '.png'))
    # plt.show()
    plt.close()

    print("E_diff > epsilon :",
          np.where(E_diff > epsilon)[0][:5], "\nepsilon: ", epsilon)
    print("ir_diff: ", np.where(ir_diff > 0)[0][:5])
    # Cache the result
    np.save('cache', np.where(E_diff > epsilon)[0][1])
    return np.where(E_diff > epsilon)[0][1]
