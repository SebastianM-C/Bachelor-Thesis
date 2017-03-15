#!/usr/bin/env python

import os
import re
import string
import numpy as np


def format_float(n):
    '''Converts the number to int if possible'''
    return ('{:n}' if n == int(n) else '{:.10g}').format(n)


def eigv_string(eigv, color='', k=-1, sep=',\\ '):
    '''Build the eigenvector string
    Add \color{} to the eigenvector (optional)
    Add \mathbf{} to the k-th element (optional)'''

    v = ('\t\t{\\color{' if color else '') + color + ('}' if color else '') + \
        '\\begin{pmatrix}\n\t\t\t'
    v += ' \\\\ '.join((('\\mathbf{' if k == j else '') + format_float(eigv[j])
                        + ('}' if k == j else ''))
                       for j in range(eigv.size))
    v += '\n\t\t\\end{pmatrix}' + ('}' if color else '') + sep + '\n'
    return v


def np_eig():
    '''Compute eigenvalues and eigenvectors with NumPy'''
    eigval, eigvec = np.linalg.eigh(H)
    eigvec = np.transpose(eigvec)
    # If the colormap element corresponding to the i-th eigenvalue is empty,
    # skip adding \color
    print(', '.join(('\\color{' + colormap[i] + '}' if colormap[i] else '') +
                    format_float(eigval[i]) for i in range(eigval.size)))
    print('\n')
    for i in range(int(nn)):
        print("v_{" + str(int(i + 1)) + '} =')
        print(eigv_string(eigvec[i], color=colormap[i],
                          k=np.argmax(np.abs(eigvec[i]))))


# Hamiltonian parameters
A = 1
B = [0, 0.1 * A, 0.2 * A, 0.4 * A, 0.5 * A]
D = [0, 0.4 * A, 0.5 * A, 0.6 * A]
N = [2, 3, 4, 5, 10, 15]

# Input parameters
b = 0
d = 0
n = 2
nn = N[n] * (N[n] + 1) / 2

# All available colors
colors = ('black', 'red', 'teal', 'blue', 'orange', 'olive',
          'magenta', 'cyan', 'Brown', 'Goldenrod', 'Green', 'Violet')

colormap = [''] * int(nn)   # colors used

# skip_zero = True
try:
    os.chdir("../Output/B" + str(B[b]) + " D" + str(D[d]) + " N" + str(N[n]))
except FileNotFoundError as fnf:
    try:
        print(fnf.strerror + ':' + fnf.filename + '\nTrying again')
        os.chdir("../../Scripts")
        os.chdir("../Output/B" + str(B[b]) +
                 " D" + str(D[d]) + " N" + str(N[n]))
        print('Succes\n' + 'Now in: ' + os.getcwd())
    except FileNotFoundError:
        print('Check the filename\n' + 'Now in: ' + os.getcwd())

# Load files
eigenvectors = np.loadtxt("eigenvectors.out", unpack=True)
H = np.loadtxt("hamilt.out")    # Hamiltonian
E = np.loadtxt("hamilt.dat", usecols=1, unpack=True)    # energy levels
index = np.loadtxt("index.out")

# Transpose matrices
H = np.transpose(H)
eigenvectors = np.transpose(eigenvectors)   # each eigenvector is on one row
eigenvalues = ''
# Group energy levels such that a level contains all the eigenvalues with
# the same value
levels = np.split(E, np.where(np.diff(E) != 0)[0] + 1)
k = 0
for i in range(len(levels)):
    for j in range(levels[i].size):
        if i < len(colors):
            colormap[i + j + k] = colors[i]
        else:
            colormap[i + j + k] = ''
    k += levels[i].size - 1
# Build eigenvalue string
# If the colormap element corresponding to the i-th eigenvalue is empty,
# skip adding \color
eigenvalues = ', '.join(
    ('\\color{' + colormap[i] + '}' if colormap[i] else '') +
    format_float(E[i]) for i in range(E.size))

# np_eig()

with open("eigenvectors.tex", "w") as f:
    # Set paper size according to the imput
    if n < 4:
        paper_size = str(6 - n)
        if n < 2:
            paper_size = str(4)
        f.write("\\documentclass[a" + paper_size +
                "paper,12pt,landscape]{article}\n\n")
        # Number of eigenvectors on a single row
        no_eigv = 10 if b == 0 and d == 0 else 5
    else:
        f.write("\\documentclass[12pt,landscape]{article}\n")
        f.write("\\setlength{\paperwidth}{" + str(int(2**(n - 3))) + "00cm}\n")
        f.write("\\setlength{\paperheight}{" + str(int(2**(n - 4))) + "00cm}")
        no_eigv = 20
    # LaTeX packages
    f.write("\\usepackage[margin=1cm]{geometry}\n")
    f.write("\\usepackage{amsmath,amsfonts,amssymb}\n")
    f.write("\\usepackage{physics}")
    f.write("\\usepackage[dvipsnames]{xcolor}\n")
    # Additional options
    f.write("\\setcounter{MaxMatrixCols}{200}\n")
    f.write("\\allowdisplaybreaks\n\n")
    # The document begins here
    f.write("\\begin{document}\n\n")

    f.write("\tEigenvectors:\n")
    f.write("\t\\begin{align*}\n")
    for i in range(int(nn)):
        f.write(('\\\\' if i % no_eigv == 0 and i != 0 else '') +
                '\t\tv_{' + str(i + 1) + '} ' +
                ('&' if i % no_eigv == 0 else '') +
                '=\n')
        # The index of the largest coefficient
        c_max = np.argmax(eigenvectors[i])
        # The energy corresponding to the eigenvector is given by the
        # (i+1)-th eigenvalue
        # Test if the eigenvectors correspond to the eigenvalues
        if b == 0 and d == 0 and \
                np.where(np.abs(H[c_max] - E[i]) < 1e-6)[0].size == 0:
            print('Warning: v' + str(i + 1) +
                  ' does not correspond to the eigenvalue ' + str(E[i]))

        if i != int(nn) - 1:
            f.write(eigv_string(eigenvectors[i], colormap[i], c_max))
        else:   # now new line after the last eigenvector
            f.write(eigv_string(eigenvectors[i], colormap[i], c_max, sep=''))
    f.write("\t\\end{align*}\n")

    f.write("\tHamiltonian:\n")
    f.write("\t\\[\n\tH=\n\t\\begin{pmatrix}\n")
    for line in H:
        f.write("\t\t")
        f.write(' & '.join(format_float(line[i]) for i in range(line.size)))
        f.write('\\\\\n')
    f.write("\t\\end{pmatrix}\n\t\\]\n")

    f.write("\tEigenvalues:\n")
    f.write("\t\\[\lambda = " + eigenvalues + '\\]\n')

    f.write("\tStates: \\[")
    for i in index:
        f.write("\\ket{" + format_float(i[0]) +
                '\, ' + format_float(i[1]) + "}")
    f.write("\\]")
    # End of the document
    f.write("\n\n\\end{document}")

os.chdir("../../Scripts")
