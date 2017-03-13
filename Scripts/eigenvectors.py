#!/usr/bin/env python

import os
import re
import string
import numpy as np


def format_float(n):
    return ('{:n}' if n == int(n) else '{:.10g}').format(n)


def eig(i, color):
    if not color:
        v = "\t\t\\begin{pmatrix}\n\t\t\t"
        v += ' \\\\ '.join(format_float(j) for j in eigenvectors[i])
        v += "\n\t\t\\end{pmatrix},\\ \n"
    else:
        v = '\t\t{\\color{' + color + '}\\begin{pmatrix}\n\t\t\t'
        v += ' \\\\ '.join(format_float(j) for j in eigenvectors[i])
        v += "\n\t\t\\end{pmatrix}},\\ \n"
    return v


# Hamiltonian parameters
A = 1
B = [0, 0.1 * A, 0.2 * A, 0.4 * A, 0.5 * A]
D = [0, 0.4 * A, 0.5 * A, 0.6 * A]
N = [2, 3, 4, 5, 10, 15]

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
index = np.loadtxt("index.out")
eigenvalues = ''
# group energy levels such that a level contains all the eigenvalues with
# the same value
levels = np.split(E, np.where(np.diff(E) != 0)[0] + 1)
k = 0
for i in range(len(levels)):
    for j in range(levels[i].size):
        if i < len(colors):
            colormap[i + j + k] = colors[i]
        else:
            colormap[i + j + k] = colors[0]
    k += levels[i].size - 1

eigenvalues = ', '.join(
    '\\color{' + colormap[i] + '}' + format_float(E[i]) for i in range(E.size))

with open("eigenvectors.tex", "w") as f:
    if n < 4:
        f.write("\\documentclass[a2paper,12pt,landscape]{article}\n\n")
    else:
        f.write("\\documentclass[12pt,landscape]{article}\n")
        f.write("\\setlength{\paperwidth}{" + str(int(2**(n - 3))) + "00cm}\n")
        f.write("\\setlength{\paperheight}{" + str(int(2**(n - 4))) + "00cm}")
    f.write("\\usepackage[margin=1cm]{geometry}\n")
    f.write("\\usepackage{amsmath,amsfonts,amssymb}\n")
    f.write("\\usepackage{physics}")
    f.write("\\usepackage[dvipsnames]{xcolor}\n")
    f.write("\\setcounter{MaxMatrixCols}{200}\n")
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
    f.write(' \\\\ '.join(format_float(j) for j in eigenvectors[i]))
    f.write("\n\t\t\\end{pmatrix}} \n")
    f.write("\t\\end{align*}\n")

    f.write("\tHamiltonian:\n")
    f.write("\t\\[\nH=\n\\begin{pmatrix}\n")
    for line in H:
        f.write("\t\t")
        f.write(' & '.join(format_float(line[i]) for i in range(line.size)))
        f.write('\\\\\n')
    f.write("\t\\end{pmatrix}\n\\]\n")

    f.write("\tEigenvalues:\n")
    f.write("\t\\[\lambda = " + eigenvalues + '\\]\n')

    f.write("\tStates: \\[")
    for i in index:
        f.write("\\ket{" + format_float(i[0]) +
                '\, ' + format_float(i[1]) + "}")
    f.write("\\]")

    f.write("\n\n\\end{document}")

os.chdir("../../Scripts")
