#!/usr/bin/env python

import os
import re
import string
import numpy as np


def eig(i):
    v = "\t\t\\begin{pmatrix}\n\t\t\t"
    v += re.sub('[\[\]]', '',
                np.array2string(eigenvectors[i], separator=' \\\\ '))
    v += "\n\t\t\\end{pmatrix},\\ \n"
    return v


# Hamiltonian parameters
A = 1
B = [0, 0.1 * A, 0.4 * A, 0.5 * A]
D = [0, 0.4 * A, 0.5 * A, 0.6 * A]
N = [2, 3, 4, 40, 60, 80, 100, 120, 140, 160]

# Input parameters
b = 0
d = 0
n = 2
nn = N[n] * (N[n] + 1) / 2

os.chdir("../Output/B" + str(B[b]) + " D" + str(D[d]) + " N" + str(N[n]))

eigenvectors = np.loadtxt("hamilt.out", unpack=True)
E = np.loadtxt("hamilt.dat", usecols=(1,), unpack=True)
eigenvalues = re.sub('[\[\]]', '', np.array2string(E, separator=', '))

with open("eigenvectors.tex", "w") as f:
    f.write("\\documentclass[a4paper,12pt]{article}\n\n")
    f.write("\\usepackage{amsmath,amsfonts,amssymb}\n\n")
    f.write("\\begin{document}\n\n")
    f.write("\tEigenvectors:\n")
    f.write("\t\\begin{multline*}\n")
    f.write("\t\tv_1" + ' =\n')
    f.write(eig(0))
    for i in range(int(nn - 2)):
        if i % 5 == 0 and i != 0:
            f.write("\\\\")
        f.write("\t\tv_" + str(i + 2) + ' =\n')
        f.write(eig(i + 1))
    # write the last eigenvector separately
    f.write("\t\tv_{" + str(int(nn)) + '} =\n')
    f.write("\t\t\\begin{pmatrix}\n\t\t\t")
    f.write(re.sub('[\[\]]', '',
                   np.array2string(eigenvectors[int(nn) - 1],
                                   separator=' \\\\ ')))
    f.write("\n\t\t\\end{pmatrix} \n")
    f.write("\t\\end{multline*}\n")
    f.write("\tEigenvalues:\n")
    f.write("\t\\[\lambda = " + eigenvalues + '\\]')
    f.write("\n\n\\end{document}")

os.chdir("../../Scripts")
