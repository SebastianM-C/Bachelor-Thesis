#!/usr/bin/env python

import os
import numpy as np
from timeit import default_timer as timer

import eigensystem
import diff
from custom_parser import parse
from tools import cd


def main(b, d, n, use_sc, stable_only=True):
    print("Running with: B = " + str(b) + " D = " + str(d) + " N = " + str(n))
    cd(b, d, n)
    start = timer()
    E, ket = eigensystem.get(use_sc)
    ir_reps = eigensystem.levels(E, ket, use_sc)

    # Select only the stable levels
    if stable_only:
        stable_levels = diff.stable(E, ir_reps, b, d, n, use_sc)
        E = E[:stable_levels]
    stop = timer()
    print("get data: ", stop - start)
    rebde = open('rebde2' + ('_sc.dat' if use_sc else '.dat'), 'w')
    reuna = open('reuna2' + ('_sc.dat' if use_sc else '.dat'), 'w')
    reuns = open('reuns2' + ('_sc.dat' if use_sc else '.dat'), 'w')

    # Write only one level corresponding to the bidimensional representation
    skip_next = False
    for i in range(E.size):
        n1 = ket[i][0]
        n2 = ket[i][1]
        if ir_reps[i] == 2:     # bidimensional representation (rebde)
            if not skip_next:
                rebde.write(str(E[i]) + '\t' + str(n1) + '\t' + str(n2) + '\n')
                skip_next = True
            else:
                skip_next = False
        else:
            if n2 % 2:          # unidimensional anti-symmetric representation
                reuna.write(str(E[i]) + '\t' + str(n1) + '\t' + str(n2) + '\n')
            else:               # unidimensional symmetric representation
                reuns.write(str(E[i]) + '\t' + str(n1) + '\t' + str(n2) + '\n')

    rebde.close()
    reuna.close()
    reuns.close()

    os.chdir("../../Scripts")
    print("Done")


if __name__ == '__main__':
    B, D, N, use_sc = parse()
    for b in B:
        for d in D:
            for n in N:
                main(b, d, n, use_sc)
