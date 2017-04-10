#!/usr/bin/env python

import os
import numpy as np
from timeit import default_timer as timer

import tools
import diff
import eigensystem


def color(x, exists):
    if x == 1 and exists:
        return "\\cellcolor{yellow}"
    if x == 2 and exists:
        return "\\cellcolor{blue!50}"
    if x == 3 and exists:
        return "\\cellcolor{green!80}"
    if not exists:
        return "\\cellcolor{red}"
    return ''


# Hamiltonian parameters
# A = 1
# B = [0, 0.1 * A, 0.2 * A, 0.4 * A, 0.5 * A]
# D = [0, 0.4 * A, 0.5 * A, 0.6 * A]
# N = [40, 60, 80, 100, 120, 140, 160]

# Input parameters
b = 0.2
d = 0.4
n = 120
stable_only = True   # choose all levels or only the stable ones
use_sc = True

start = timer()
# os.chdir("../Output/B" + str(B[b]) + " D" + str(D[d]) + " N" + str(N[n]))
tools.cd(b, d, n)

# Read data
# if stable_only:
#     E = np.loadtxt("s_hamilt.dat", usecols=(0,), unpack=True)
E, ket = eigensystem.get(use_sc)
ir_reps = eigensystem.levels(E, ket, use_sc)
if stable_only:
    delta_n = 10 if not use_sc else 20
    stable_levels = diff.stable(E, ir_reps, b, d, n, use_sc, delta_n)
    E = E[:stable_levels]

rebde = np.loadtxt("rebde2" + ('_sc.dat' if use_sc else '.dat'),
                   usecols=0, unpack=True)
reuna = np.loadtxt("reuna2" + ('_sc.dat' if use_sc else '.dat'),
                   usecols=0, unpack=True)
reuns = np.loadtxt("reuns2" + ('_sc.dat' if use_sc else '.dat'),
                   usecols=0, unpack=True)

# Bi-directional search
E_in_rebde = np.in1d(E, rebde, assume_unique=False)
E_in_reuna = np.in1d(E, reuna, assume_unique=True)
E_in_reuns = np.in1d(E, reuns, assume_unique=True)

rebde_in_E = np.in1d(rebde, E, assume_unique=False)
reuna_in_E = np.in1d(reuna, E, assume_unique=True)
reuns_in_E = np.in1d(reuns, E, assume_unique=True)

if E.size < max(rebde.size, reuna.size, reuns.size):
    rebde = rebde[:E.size]
    reuna = reuna[:E.size]
    reuns = reuns[:E.size]

    rebde_in_E = rebde_in_E[:E.size]
    reuna_in_E = reuna_in_E[:E.size]
    reuns_in_E = reuns_in_E[:E.size]

# Padding
rebde = np.pad(rebde, pad_width=(0, E.size - rebde.size), mode='constant')
reuna = np.pad(reuna, pad_width=(0, E.size - reuna.size), mode='constant')
reuns = np.pad(reuns, pad_width=(0, E.size - reuns.size), mode='constant')

rebde_in_E = np.pad(rebde_in_E, pad_width=(
    False, E.size - rebde_in_E.size), mode='constant')
reuna_in_E = np.pad(reuna_in_E, pad_width=(
    False, E.size - reuna_in_E.size), mode='constant')
reuns_in_E = np.pad(reuns_in_E, pad_width=(
    False, E.size - reuns_in_E.size), mode='constant')

files = np.array([E, rebde, reuna, reuns, E_in_rebde, E_in_reuna,
                  E_in_reuns, rebde_in_E, reuna_in_E, reuns_in_E])

with open("table.tex", "w") as f:
    f.write("\\documentclass{article}\n\n")
    f.write("\\usepackage[margin=0.2in]{geometry}\n\\usepackage{longtable}\n")
    f.write("\\usepackage[table]{xcolor}\n\n")
    f.write("\\begin{document}\n\n")
    f.write("\\begin{longtable}{" + " | ".join(["c"] * 4) + "}\n")
    f.write("hamilt.dat & rebde.dat & reuna.dat & reuns.dat\t\\\\\n")
    f.write("\\hline\n\\endfirsthead\n")
    for row in range(files.shape[1]):
        # Find the representation of the energy level
        c = 0
        for x in range(1, 4):
            if files[x + 3][row]:
                c = x

        line = color(c, True) + str(files[0][row]) + " & " + \
            " & ".join(color(x, files[x + 6][row]) +
                       str(files[x][row]) for x in range(1, 4))
        f.write('\t' + line + " \\\\\n")
    f.write("\\end{longtable}")
    f.write("\n\n\\end{document}")

os.chdir("../../Scripts")
end = timer()
print(end - start)
