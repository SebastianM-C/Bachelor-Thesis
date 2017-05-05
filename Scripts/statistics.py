#!/usr/bin/env python

import os
import numpy as np

import select_rep
from tools import cd
from custom_parser import parse
from diff import relSpacing
from plots import bar_plot, histogram


def main(b, d, n, delta_n, st_epsilon, lvl_epsilon, reselect=True, cut=0,
         bin_size=0.25):
    if reselect:
        select_rep.main(b, d, n, delta_n, st_epsilon, lvl_epsilon, cut)
    reps = 'reuna', 'reuns', 'rebde'
    cd(b, d, n)
    deltaE = np.loadtxt('stable.txt')[1]
    count = int(4 / bin_size) + 1
    rel_sp = []
    w = []          # weights
    for r in reps:
        rep = np.loadtxt(r + '.dat', usecols=(0,))
        rel_sp.append(relSpacing(rep))
        w.append(np.ones(rel_sp[-1].shape) / 3)
        histogram(rel_sp[-1], bins=np.linspace(0, 4, count), label=r,
                  fname=r + '.png', xlabel='S')
        bar_plot(rel_sp[-1], label=r, ylabel='S',
                 fname='bar_P(S)_' + r + '.png',
                 dpi=400)

    # Relative spacing histogram
    histogram(rel_sp, bins=np.linspace(0, 4, count), weights=w,
              normed=True, label=reps, ylabel='$P(S)$', xlabel='$S$',
              fname='P(S)' + '_st_' + '{:.0e}'.format(st_epsilon) + '_eps_' +
              '{:.0e}'.format(lvl_epsilon) +
              ('_cut_' + '{:.2f}'.format(cut) + '.png' if cut else '.png'),
              count=count, stacked=True, use_wigner=True, use_poisson=True)
    histogram(rel_sp, bins=np.linspace(0, 4, count), weights=w,
              normed=True, title='$\\Delta E =' + '{:.5}'.format(deltaE) + '$',
              fname='P(S)' + '_fit_' + '{:.0e}'.format(st_epsilon) + '_eps_' +
              '{:.0e}'.format(lvl_epsilon) +
              ('_cut_' + '{:.2f}'.format(cut) + '.png' if cut else '.png'),
              count=count, ylabel='$P(S)$', xlabel='$S$', stacked=True,
              label=reps, fit=True)
    histogram(rel_sp, cumulative=True, bins=np.linspace(0, 4, count),
              normed=True, ylabel=r'$\sum_0^S P(x)$', xlabel='$S$', label=reps,
              use_wigner=True, count=count, fname='Cumulative P(S).png')

    os.chdir("../../Scripts")


if __name__ == '__main__':
    B, D, N, delta_n, st_epsilon, lvl_epsilon, reselect, cut, bin_size = \
        parse(advanced=True, select=True, hist_bin=True)

    for b in B:
        for d in D:
            for n in N:
                main(b, d, n, delta_n, st_epsilon, lvl_epsilon, reselect, cut,
                     bin_size)
