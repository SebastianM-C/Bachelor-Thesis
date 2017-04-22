#!/usr/bin/env python

import os
import numpy as np
# import matplotlib.pyplot as plt

import select_rep
from tools import cd
from custom_parser import parse
from diff import relSpacing
from plots import bar_plot, histogram


def main(b, d, n, use_sc, delta_n, st_epsilon, lvl_epsilon, reselect=True,
         cut=0):
    if reselect:
        select_rep.main(b, d, n, use_sc, delta_n, st_epsilon, lvl_epsilon, cut)
    reps = 'reuna', 'reuns', 'rebde'
    cd(b, d, n)
    rel_sp = []
    w = []          # weights
    for r in reps:
        rep = np.loadtxt(r + '2' + ('_sc.dat' if use_sc else '.dat'),
                         usecols=0)
        rel_sp.append(relSpacing(rep))
        w.append(np.ones(rel_sp[-1].shape) / 3)
        histogram(rel_sp[-1], bins=np.arange(0, 4, 1 / 4), label=r,
                  fname=r + ('_sc.png' if use_sc else '.png'), show=False)
        bar_plot(rel_sp[-1], label=r, ylabel='S',
                 fname='bar_P(S)_' + r + ('_sc.png' if use_sc else '.png'),
                 dpi=400)

    # Relative spacing histogram
    histogram(rel_sp, bins=np.arange(0, 4, 1 / 4), weights=w, stacked=True,
              normed=True, label=reps, ylabel='P(S)', xlabel='S',
              fname='P(S)' + '_st_' + '{:.0e}'.format(st_epsilon) + '_eps_' +
              '{:.0e}'.format(lvl_epsilon) +
              ('_cut_' + '{:.2f}'.format(cut) + '.png' if cut else '.png'),
              use_wigner=True, use_poisson=True)
    histogram(rel_sp, cumulative=True, bins=np.arange(0, 4, 1 / 4), label=reps,
              normed=True, ylabel='P(S)', xlabel='S', use_wigner=True,
              fname='Cumulative P(S)' + ('_sc.png' if use_sc else '.png'))

    os.chdir("../../Scripts")


if __name__ == '__main__':
    B, D, N, use_sc, delta_n, st_epsilon, lvl_epsilon, reselect, cut = \
        parse(advanced=True, select=True)

    for b in B:
        for d in D:
            for n in N:
                main(b, d, n, use_sc, delta_n, st_epsilon, lvl_epsilon,
                     reselect, cut)
