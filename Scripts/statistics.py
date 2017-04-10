#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt

import select_rep
from tools import cd
from custom_parser import parse


def relSpacing(rep):
    deltaE = np.diff(rep)
    avgSpacing = (rep[-1] - rep[0]) / rep.size
    return deltaE / avgSpacing


def histogram(data, name, use_sc, show, save=True):
    plt.hist(data, bins=np.arange(0, 4, 1 / 4), label=name)
    plt.legend()
    if save:
        plt.savefig(name + ('_sc.png' if use_sc else '.png'))
    if show:
        plt.show()
    plt.close()


def main(b, d, n, use_sc, re_select=True, show=False):
    if re_select:
        delta_n = 20    # if use_sc else 10
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
        histogram(rel_sp[-1], r, use_sc, show)

    plt.hist(rel_sp, bins=np.arange(0, 4, 1 / 4), weights=w, stacked=True,
             normed=True, label=reps)
    plt.ylabel('P(S)')
    plt.xlabel('S')
    plt.legend()
    plt.savefig('P(S)' + ('_sc.png' if use_sc else '.png'))
    # plt.show()
    plt.close()
    plt.hist(rel_sp, cumulative=True, bins=np.arange(0, 4, 1 / 4),
             normed=True, label=reps)
    plt.ylabel('P(S)')
    plt.xlabel('S')
    plt.legend()
    plt.savefig('Cumulative P(S)' + ('_sc.png' if use_sc else '.png'))
    # plt.show()
    plt.close()

    os.chdir("../../Scripts")


if __name__ == '__main__':
    B, D, N, use_sc = parse()
    for b in B:
        for d in D:
            for n in N:
                main(b, d, n, use_sc)
