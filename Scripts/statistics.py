#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt

import select_rep
from tools import cd
from custom_parser import parse
from diff import relSpacing
from eigensystem import get_state
from plots import bar_plot, histogram


# b, d, n = 0.2, 0.4, 60
def main(b, d, n, use_sc, re_select=True, show=False):
    if re_select:
        delta_n = 20
        select_rep.main(b, d, n, use_sc, delta_n)
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
                  fname=r + ('_sc.png' if use_sc else '.png'), show=show)
        bar_plot(rel_sp[-1], label=r, ylabel='S',
                 fname='bar_P(S)_' + r + ('_sc.png' if use_sc else '.png'),
                 dpi=400)

    # Relative spacing histogram
    histogram(rel_sp, bins=np.arange(0, 4, 1 / 4), weights=w, stacked=True,
              normed=True, label=reps, ylabel='P(S)', xlabel='S',
              fname='P(S)' + ('_sc.png' if use_sc else '.png'),
              use_wigner=True)
    histogram(rel_sp, cumulative=True, bins=np.arange(0, 4, 1 / 4), label=reps,
              normed=True, ylabel='P(S)', xlabel='S',
              fname='Cumulative P(S)' + ('_sc.png' if use_sc else '.png'))
    # # Energy difference (between two consecutive levels) histogram
    # state, _, _ = get_state(use_sc)
    # stable = int(np.loadtxt('cache.txt'))
    # delta = np.diff(state['E'])[:stable]
    # epsilon = 1e-5
    # histogram(delta, label='$\\Delta = E_{n+1} - E_n$ stable',
    #           bins=np.pad(np.geomspace(1e-9, 1e3, 13), (1, 0),
    #                       mode='constant'), xscale='log',
    #           fname='hist_delta_stable' + ('_sc.png' if use_sc else '.png'))
    # # Energy difference bar plot
    # bar_plot(delta, label='$\\Delta = E_{n+1} - E_n$ stable',
    #          figsize=(20, 4), yscale='log', axhline_y=epsilon, dpi=600,
    #          fname='bar_delta_stable' + ('_sc.png' if use_sc else '.png'),
    #          bbox_inches='tight')

    os.chdir("../../Scripts")


if __name__ == '__main__':
    B, D, N, use_sc = parse()
    for b in B:
        for d in D:
            for n in N:
                main(b, d, n, use_sc)
