#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt


def ir_decomp(E):
    "Irreductible representations decomposition"
    reuna = []
    reuns = []
    reuns.append(E[0])

    # in order to compute E[j] - E[j-1] using np.diff() we must shift every
    # element of the array to right with 1 position and ignore diff[0]
    diff = np.diff(np.roll(E, 1))           # E[j] - E[j-1]
    E = np.roll(E, 1)
    E = E[1:]
    diff1 = np.diff(E)                      # E[j+1] - E[j]

    # print("size of E: ", E.size, "\nsize of diff: ", diff.size)
    print("diff: ", diff[:20])

    rebde = E[diff <= 0.0005]
    rebde = rebde[1:]
    print("\nrebde: ", rebde[:10], "...", rebde[-10:], rebde.size)

    return (rebde, reuna, reuns)


def statistics(A, B, D, N, f_index):
    "Compute spacing statistics"
    path = "../Output/B" + str(B) + " N" + str(N)
    os.chdir(path)
    output = open("diff.dat", 'w')

    E = np.loadtxt(input[0], usecols=(1,), unpack=True)  # Energy levels

    # Select only stable levels
    E = E[stable_levels]

    if f_index > 0 & f_index < 4:
        # The energy levels must be decomposed again in
        # irreductible representations
        E = ir_decomp(E)[f_index - 1]
    #    print(ir_rep)

    deltaE = np.diff(E)                                  # Energy level spacing
    avgSpacing = (E[-1] - E[0]) / E.size                 # Average spacing
    print("\nStatistics parameters: \nA = ", A, "\nB = ", B, "\nD = ",
          D, "\nN = ", N, "\nAverage spacing: ", avgSpacing)
    relSpacing = deltaE / avgSpacing
    relSpacing.tofile(output, sep="\n")
    plt.hist(relSpacing, bins=np.arange(0, 4, 1 / 2))
    plt.show()
    os.chdir("../../Scripts")
    return (deltaE, avgSpacing, relSpacing)


def difference(A, B, D, N1, N2, file):
    "Compare the energy levels between different diagonalization bases"
    os.chdir("../Output")
    path1 = "B" + str(B) + " N" + str(N1)
    path2 = "B" + str(B) + " N" + str(N2)
    os.chdir(path1)

    E1 = np.loadtxt(file, usecols=(1,), unpack=True)  # Energy levels
    os.chdir("../" + path2)
    E2 = np.loadtxt(file, usecols=(1,), unpack=True)  # Energy levels

    if E2.size > E1.size:
        diff = (E2[:E1.size] - E1) / E1
    else:
        diff = (E1[:E2.size] - E2) / E2

    os.chdir("../../Statistics")
    output = open("E_diff.dat", 'w')
    diff.tofile(output, sep="\n")
    os.chdir("../Scripts")
    return diff


# Hamiltonian parameters
A = 1
B = [0, 0.1 * A, 0.5 * A]
D = [0, 0.4 * A, 0.5 * A]
N = [40, 60, 80, 100, 120, 140, 160]
input = ["hamilt.dat", "rebde.dat", "reuna.dat", "reuns.dat"]

# Input parameters
b = 0
d = 2
i = 1       # first diagonalization basis size
j = 2       # second diagonalization basis size

# Enargy level difference at a change of basis
diff = difference(A, B[b], D[d], N[i], N[j], input[0])
stable_levels = np.where(np.abs(diff) < 5e-2)
print(diff[np.abs(diff) < 5e-2].size,
      " energy levels varry with less than 5% in a change of basis from N =",
      N[i], "to N =", N[j])
print("Stable levels: ", stable_levels)
# only 0 and 1 are valid for the last argument; can't separate reuna and reuns
stats = statistics(A, B[b], D[d], N[i], 1)

print("\nEnargy level difference at a change of basis: \n")
plt.hist(diff, bins='auto')
plt.show()
