#!/usr/bin/env python

import os
import re
import string
import numpy as np
from scipy import linalg
from timeit import default_timer as timer


def format_float(n):
    '''Converts the number to int if possible'''
    return ('{:n}' if n == int(n) else '{:.8g}').format(n)


def optional_color(element, color):
    '''Add optional color to the element'''
    return '{\\color{' + color + '}' + element + '}' if color else element


def eigv_string(eigv, color='', k=-1, sep=',\\ '):
    '''Build the eigenvector string
    Add \color{} to the eigenvector (optional)
    Add \mathbf{} to the k-th element (optional)'''

    v = ('\t\t{\\color{' + color + '}' if color else '') + \
        '\\begin{pmatrix}\n\t\t\t'
    # Limit the number of displayed eigenvector elements
    eigv_size = eigv.size if eigv.size < max_display else int(max_display / 4)
    v += ' \\\\ '.join(('\\mathbf{' + format_float(eigv[j]) + '}' if k == j
                        else format_float(eigv[j]))
                       for j in range(eigv_size))
    v += '\n\t\t\\end{pmatrix}' + ('}' if color else '') + sep + '\n'
    return v


def h_elem(elem, color='', skip_zero=False):
    '''Build Hamiltonian element string
    Optionally add color and skip zeros'''
    return (optional_color(format_float(elem), color)) \
        if (not skip_zero) or elem != 0 else ''


def ket(n1, n2, color='', bf=False):
    '''Build \ket{} with optional color and optional bold face'''
    k = '\\ket{' + str(n1) + '\, ' + str(n2) + '}'
    k = '\\mathbf{' + k + '}' if bf else k
    return optional_color(k, color)


# Hamiltonian parameters
A = 1
B = [0, 0.1 * A, 0.2 * A, 0.4 * A, 0.5 * A, 0.6 * A, 0.8 * A]
D = [0, 0.4 * A, 0.5 * A, 0.6 * A]
N = [2, 3, 4, 5, 6, 10, 15, 50]

# Input parameters
b = 2
d = 1
n = 2
nn = int(N[n] * (N[n] + 1) / 2)

# Use hamilt or SciPy for eigenvalues and eigenvectors
use_sc = False

# Maximum number of elements to display
max_display = 300

# All available colors
colors = ('black', 'red', 'teal', 'blue', 'orange', 'olive',
          'magenta', 'cyan', 'Brown', 'Goldenrod', 'Green', 'Violet')

colormap = [''] * nn   # colors used

dgc = b == 0 and d == 0     # degenerate case

# Try to cd to the given path. In case of an error go back to ../../Scripts
# and try again (maybe the last run had an error or
# the script did not reach the end)
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
H = np.loadtxt("hamilt.out")    # Hamiltonian
index = np.loadtxt("index.out", dtype=int)
c_max = np.empty([nn], dtype=int)

# Transpose matrices
H = np.transpose(H)
# Get eigenvalues and eigenvectors
if use_sc:
    E, eigenvectors = linalg.eigh(H)
else:
    E = np.loadtxt("hamilt.dat", usecols=1, unpack=True)    # energy levels
    eigenvectors = np.loadtxt("eigenvectors.out", unpack=True)

eigenvectors = np.transpose(eigenvectors)   # each eigenvector is on one row
eigenvalues = ''
# Group energy levels such that a level contains all the eigenvalues with
# the same value
epsilon = 5e-4 if not use_sc else 1e-6
levels = np.split(E, np.where(np.diff(E) > epsilon)[0] + 1)

k = 0
for i in range(len(levels)):
    for j in range(levels[i].size):
        if i < len(colors):
            colormap[i + j + k] = colors[i]
        else:
            colormap[i + j + k] = ''
    k += levels[i].size - 1

# The index of the largest coefficient
for i in range(nn):
    c_max[i] = np.argmax(np.abs(eigenvectors[i]))
# Build eigenvalue string
# If the colormap element corresponding to the i-th eigenvalue is empty,
# skip adding \color
# Limit the number of eigenvalues displayed
E_size = E.size if E.size < max_display else int(max_display / 2)
eigenvalues = ', '.join(
    optional_color(format_float(E[i]), colormap[i]) +
    ('\n\t' if i % 5 == 0 and i else '')   # add newline for shorter lines
    for i in range(E_size))

with open("results B" +
          str(B[b]) + ' D' + str(D[d]) + ' N' + str(N[n]) +
          ('_sc' if use_sc else '') +
          ".tex", "w") as f:
    # Set paper size according to the imput
    if n < 5:
        # Set paper size and number of eigenvectors on a single row
        if n < 3:
            paper_size = str(4)
            no_eigv = 5
        if n == 3:
            paper_size = str(2)
            no_eigv = 10
        if n == 4:
            paper_size = str(1)
            no_eigv = 15
        if dgc:
            no_eigv *= 2

        f.write("\\documentclass[a" + paper_size +
                "paper,12pt,landscape]{article}\n\n")
    else:
        paper_w = int(2**(n - 3)) * 50
        paper_h = int(2**(n - 4)) * 50
        # Limit paper size
        if paper_w > 500:
            paper_w = 350
        if paper_h > 500:
            paper_h = 300
        f.write("\\documentclass[12pt,landscape]{article}\n")
        f.write("\\setlength{\paperwidth}{" + str(paper_w) + 'cm' + "}\n")
        f.write("\\setlength{\paperheight}{" + str(paper_h) + 'cm' + "}\n")
        no_eigv = (n - 3) * 20 + 10 \
            if nn < max_display else int(max_display / 5)

    # LaTeX packages
    f.write("\\usepackage[margin=1cm]{geometry}\n")
    f.write("\\usepackage{amsmath,amsfonts,amssymb}\n")
    f.write("\\usepackage{physics}")
    f.write("\\usepackage[dvipsnames]{xcolor}\n")
    # Additional options
    f.write("\\setcounter{MaxMatrixCols}{" +
            (str(nn + 10) if nn < max_display else str(max_display)) + "}\n")
    f.write("\\allowdisplaybreaks\n\n")
    # The document begins here
    f.write("\\begin{document}\n\n")

    f.write("\tEigenvectors:\n")
    f.write("\t\\begin{align*}\n")

    for i in range(nn):
        # Limit the number of displayed eigenvectors
        if i == max_display:
            f.write("\t\t\cdots")
            break
        f.write(('\\\\' if i % no_eigv == 0 and i != 0 else '') +
                '\t\tv_{' + str(i + 1) + '} ' +
                ('&' if i % no_eigv == 0 else '') +
                '=\n')

        # The energy corresponding to the eigenvector is given by the
        # (i+1)-th eigenvalue
        # Test if the eigenvectors correspond to the eigenvalues
        # in the degenerate case
        if dgc and np.where(np.abs(H[c_max[i]] - E[i]) < 1e-6)[0].size == 0:
            print('Warning: v' + str(i + 1) +
                  ' does not correspond to the eigenvalue ' + str(E[i]))

        if i != nn - 1:
            f.write(eigv_string(eigenvectors[i], colormap[i], c_max[i]))
        else:   # now new line after the last eigenvector
            f.write(eigv_string(eigenvectors[i],
                                colormap[i], c_max[i], sep=''))

    f.write("\t\\end{align*}\n")

    # Skip displaying the Hamiltonian if the size is too large
    if nn < max_display:
        f.write("\tHamiltonian:\n")
        f.write("\t\\[\n\tH=\n\t\\begin{pmatrix}\n")
        for i in range(nn):
            f.write("\t\t")
            # Add colors on the diagonal and remove zeros in the degenerate
            # case
            f.write(' & '.join(
                h_elem(H[i][j]) if not dgc else
                h_elem(H[i][j],
                       color=colors[int(round(H[i][j]))] if i == j else '',
                       skip_zero=True if not (i == j) else False)
                    for j in range(H[i].size)))
            f.write('\\\\\n')
        f.write("\t\\end{pmatrix}\n\t\\]\n")

    f.write("\tEigenvalues:\n")
    f.write("\t\\[\lambda = " + eigenvalues + '\\]\n')

    f.write("\tStates: \n\t\\[")
    # Limit the number of displayed elements
    for i in index[:max_display]:
        f.write(ket(i[0], i[1]))
    f.write("\\]\n")

    f.write("\tOrdering: \n\t\\[")
    for i in range(nn):
        # Limit the number of displayed elements
        if i == max_display:
            f.write("\t\t\cdots")
            break
        n1 = index[c_max[i]][0]
        n2 = index[c_max[i]][1]
        f.write(ket(n1, n2, colormap[i]) + '\, ' +
                # add newline for shorter lines
                ('\n\t' if i % 4 == 0 and i else '')
                )
    f.write("\\]\n")

    # Check for duplicate states
    unique, counts = np.unique(c_max, return_counts=True)
    if np.any(counts > 1):
        # print("Warning: Duplicates")
        f.write("\tDuplicates: \n\t\\[")
        for i in (index[unique[counts > 1]])[:max_display]:
            f.write(ket(i[0], i[1]))
        f.write("\\]\n")

    # End of the document
    f.write("\n\n\\end{document}")

os.chdir("../../Scripts")
