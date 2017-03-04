#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt


def statistics(A, B, D, N, file):
    "Compute spacing statistics"
    path = "../Output/B" + str(B) + " D" + str(D) + " N" + str(N)
    os.chdir(path)
    output = open("diff.dat", 'w')
    stable = open("s_hamilt.dat", "w")

    E = np.loadtxt(input[0], usecols=(1,), unpack=True)  # Energy levels
    # Read the chosen representation
    rep = np.loadtxt(file, usecols=(1,), unpack=True)

    # Select only stable levels
    E = E[stable_levels]
    E.tofile(stable, sep='\n')
    # Select only the stable levels in the representation
    rep = np.intersect1d(E, rep, assume_unique=True)

    deltaE = np.diff(rep)                                # Energy level spacing
    avgSpacing = (rep[-1] - rep[0]) / rep.size           # Average spacing
    print("\nStatistics parameters: \nfile: ", file, "\nA = ", A, "\nB = ", B,
          "\nD = ", D, "\nN = ", N, "\nAverage spacing: ", avgSpacing)
    relSpacing = deltaE / avgSpacing
    relSpacing.tofile(output, sep="\n")
    plt.hist(relSpacing, bins=np.arange(0, 4, 1 / 4))
    plt.show()
    os.chdir("../../Scripts")
    return (deltaE, avgSpacing, relSpacing)


def difference(A, B, D, N1, N2, file):
    "Compare the energy levels between different diagonalization bases"
    os.chdir("../Output")
    path1 = "B" + str(B) + " D" + str(D) + " N" + str(N1)
    path2 = "B" + str(B) + " D" + str(D) + " N" + str(N2)
    os.chdir(path1)

    E1 = np.loadtxt(file, usecols=(1,), unpack=True)  # Energy levels
    os.chdir("../" + path2)
    E2 = np.loadtxt(file, usecols=(1,), unpack=True)  # Energy levels

    if E2.size > E1.size:
        diff = (E2[:E1.size] - E1) / E1
    else:
        diff = (E1[:E2.size] - E2) / E2

    # os.chdir("../../Statistics")
    # output = open("E_diff.dat", 'w')
    # diff.tofile(output, sep="\n")
    os.chdir("../../Scripts")
    return diff


# Hamiltonian parameters
A = 1
B = [0, 0.1 * A, 0.4 * A, 0.5 * A]
D = [0, 0.4 * A, 0.5 * A, 0.6 * A]
N = [40, 60, 80, 100, 120, 140, 160]
input = ["hamilt.dat", "rebde.dat", "reuna.dat", "reuns.dat"]

# Input parameters
b = 3
d = 3
i = 1       # first diagonalization basis size
j = 2       # second diagonalization basis size

# Energy level difference at a change of basis
diff = difference(A, B[b], D[d], N[i], N[j], input[0])
stable_levels = np.where(np.abs(diff) < 5e-2)
print(diff[np.abs(diff) < 5e-2].size,
      " energy levels varry with less than 5% in a change of basis from N =",
      N[i], "to N =", N[j])
print("Stable levels: ", stable_levels)

stats = statistics(A, B[b], D[d], N[i], input[1])

print("\nEnergy level difference at a change of basis: \n")
plt.hist(diff, bins='auto')
plt.show()
