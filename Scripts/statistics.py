#!/usr/bin/env python

import os
import numpy as np

import select_rep
from tools import cd
from custom_parser import parse
from diff import relSpacing
from plots import bar_plot, histogram


def main(b, d, n, delta_n, st_epsilon, lvl_epsilon, reselect=True, cut=0,
         bin_size=0.25, max_energy=0):
    if reselect:
        select_rep.main(b, d, n, delta_n, st_epsilon, lvl_epsilon, cut)
    reps = 'reuna', 'reuns', 'rebde'
    cd(b, d, n)
    if max_energy:
        deltaE = max_energy
    else:
        deltaE = np.loadtxt('stable.txt')[1]
    count = int(4 / bin_size) + 1
    rel_sp = []
    avg_sp = []
    w = []          # weights
    for r in reps:
        rep = np.loadtxt(r + '.dat', usecols=(0,))
        if max_energy:
            rep = rep[rep <= max_energy]
        rel_sp.append(relSpacing(rep))
        avg_sp.append((rep[-1] - rep[0]) / rep.size)
        w.append(np.ones(rel_sp[-1].shape) / 3)
        fname = r + ('_max_e_' + str(max_energy) + '.pdf'
                     if max_energy else '.pdf')
        histogram(rel_sp[-1], bins=np.linspace(0, 4, count), label=r,
                  fname=fname, xlabel='$s$', figsize=(2.8, 3),
                  ylabel='No. of levels')
        fname = 'bar_' + r + ('_max_e_' + str(max_energy) + '.pdf'
                              if max_energy else '.pdf')
        bar_plot(rel_sp[-1], label=r, ylabel='s', xlabel='index',
                 fname=fname, dpi=400, figsize=(2.8, 3),
                 title=r'$\frac{E_n-E_0}{N}=' +
                 '{:.3}'.format(avg_sp[-1]) + '$')

    # Sace the average spacings
    fname = 'avg_sp' + \
        ('_max_e_' + str(max_energy) + '.txt' if max_energy else '.txt')
    with open(fname, 'w') as f:
        f.write('\n'.join([str(i) for i in avg_sp]))
    # Relative spacing histogram P(s)
    fname = 'P(s)' + '_st_' + '{:.0e}'.format(st_epsilon) + '_eps_' + \
        '{:.0e}'.format(lvl_epsilon) + \
        ('_cut_' + '{:.2f}'.format(cut) if cut else '') + \
        ('_max_e_' + str(max_energy) + '.pdf' if max_energy else '.pdf')
    histogram(rel_sp, bins=np.linspace(0, 4, count), weights=w,
              normed=True, label=reps, ylabel='$P(s)$', xlabel='$s$',
              fname=fname, count=count, stacked=True, use_wigner=True,
              use_poisson=True, figsize=(5.8, 4.5))
    # Fitted P(s)
    fname = 'P(s)_fit_' + '{:.0e}'.format(st_epsilon) + '_eps_' + \
        '{:.0e}'.format(lvl_epsilon) + \
        ('_cut_' + '{:.2f}'.format(cut) if cut else '') + \
        ('_max_e_' + str(max_energy) + '.pdf' if max_energy else '.pdf')
    histogram(rel_sp, bins=np.linspace(0, 4, count), weights=w,
              normed=True, title='$\\Delta E =' + '{:.5}'.format(deltaE) + '$',
              fname=fname, count=count, ylabel='$P(s)$', xlabel='$s$',
              stacked=True, label=reps, fit=True, max_e=max_energy,
              figsize=(5.8, 4.5))
    # cumulative relative spacing histogram I(s)
    fname = 'I(s)' + \
        ('_max_e_' + str(max_energy) + '.pdf' if max_energy else '.pdf')
    histogram(rel_sp, cumulative=True, bins=np.linspace(0, 4, count),
              normed=True, ylabel=r'$I(s)$', xlabel='$s$', label=reps,
              use_wigner=True, use_poisson=True, count=count, fname=fname,
              figsize=(5.8, 4.5))
    # Fitted I(s)
    fname = 'I(s)_fit' + \
        ('_max_e_' + str(max_energy) + '.pdf' if max_energy else '.pdf')
    histogram(rel_sp, cumulative=True, bins=np.linspace(0, 4, count),
              normed=True, ylabel=r'$I(s)$', xlabel='$s$', label=reps,
              fit=True, count=count, fname=fname,  figsize=(5.8, 4.5))
    # Version
    with open('version.txt', 'w') as f:
        f.write('1.3')
    os.chdir("../../Scripts")


if __name__ == '__main__':
    B, D, N, delta_n, st_epsilon, lvl_epsilon, reselect, cut, bin_size, \
        max_energy = parse(advanced=True, select=True, hist_bin=True,
                           max_e=True)

    for b in B:
        for d in D:
            for n in N:
                for max_e in max_energy:
                    main(b, d, n, delta_n, st_epsilon, lvl_epsilon, reselect,
                         cut, bin_size, max_e)
