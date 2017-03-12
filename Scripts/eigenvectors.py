#!/usr/bin/env python

import os
import re
import string
import numpy as np


def eig(i, color):
    if not color:
        v = "\t\t\\begin{pmatrix}\n\t\t\t"
        v += re.sub('[\[\]]', '',
                    np.array2string(eigenvectors[i], separator=' \\\\ '))
        v += "\n\t\t\\end{pmatrix},\\ \n"
    else:
        v = '\t\t{\\color{' + color + '}\\begin{pmatrix}\n\t\t\t'
        v += re.sub('[\[\]]', '',
                    np.array2string(eigenvectors[i], separator=' \\\\ '))
        v += "\n\t\t\\end{pmatrix}},\\ \n"
    return v


# Hamiltonian parameters
A = 1
B = [0, 0.1 * A, 0.2 * A, 0.4 * A, 0.5 * A]
D = [0, 0.4 * A, 0.5 * A, 0.6 * A]
N = [2, 3, 4, 5, 40, 60, 80, 100, 120, 140, 160]

# Input parameters
b = 0
d = 0
n = 3
nn = N[n] * (N[n] + 1) / 2

# all available colors
colors = ('black', 'red', 'teal', 'blue', 'orange', 'olive',
          'magenta', 'cyan', 'Brown', 'Goldenrod', 'Green', 'Violet')
# colors used
colormap = [''] * int(nn)

os.chdir("../Output/B" + str(B[b]) + " D" + str(D[d]) + " N" + str(N[n]))

eigenvectors = np.loadtxt("eigenvectors.out", unpack=True)
H = np.loadtxt("hamilt.out")
H = np.transpose(H)
E = np.loadtxt("hamilt.dat", usecols=(1,), unpack=True)
eigenvalues = ''
# group energy levels such that a level contains all the eigenvalues with
# the same value
levels = np.split(E, np.where(np.diff(E) != 0)[0] + 1)
k = 0
for i in range(len(levels)):
    for j in range(levels[i].size):
        colormap[i + j + k] = colors[i]
    k += levels[i].size - 1

eigenvalues = ', '.join(
    '\\color{' + colormap[i] + '}' + str(E[i]) for i in range(E.size))

with open("eigenvectors.tex", "w") as f:
    f.write("\\documentclass[a4paper,12pt,landscape]{article}\n\n")
    f.write("\\usepackage[margin=1cm]{geometry}\n")
    f.write("\\usepackage{amsmath,amsfonts,amssymb}\n")
    f.write("\\usepackage[dvipsnames]{xcolor}\n")
    f.write("\\setcounter{MaxMatrixCols}{20}")
    f.write("\\allowdisplaybreaks\n\n")
    f.write("\\begin{document}\n\n")
    f.write("\tEigenvectors:\n")
    f.write("\t\\begin{align*}\n")
    f.write("\t\tv_1" + ' &=\n')
    f.write(eig(0, ''))
    for i in range(int(nn - 2)):
        if (i + 1) % 5 == 0 and i != 0:
            f.write("\\\\")
            f.write("\t\tv_{" + str(i + 2) + '} &=\n')
        else:
            f.write("\t\tv_{" + str(i + 2) + '} =\n')
        f.write(eig(i + 1, colormap[i + 1]))
    # write the last eigenvector separately
    if nn % 5 == 0:
        f.write("\t\tv_{" + str(int(nn)) + '} =\n')
    else:
        f.write("\\\\\t\tv_{" + str(int(nn)) + '} &=\n')
    f.write(
        "\t\t{\\color{" + colormap[int(nn) - 1] + "}\\begin{pmatrix}\n\t\t\t")
    f.write(re.sub('[\[\]]', '',
                   np.array2string(eigenvectors[int(nn) - 1],
                                   separator=' \\\\ ')))
    f.write("\n\t\t\\end{pmatrix}} \n")
    f.write("\t\\end{align*}\n")
    f.write("\tHamiltonian:\n")
    f.write("\t\\[\nH=\n\\begin{pmatrix}\n")
    for line in H:
        f.write("\t\t")
        f.write(' & '.join(str(line[i]) for i in range(line.size)))
        f.write('\\\\\n')
    f.write("\t\\end{pmatrix}\n\\]\n")
    f.write("\tEigenvalues:\n")
    f.write("\t\\[\lambda = " + eigenvalues + '\\]')
    f.write("\n\n\\end{document}")

os.chdir("../../Scripts")
