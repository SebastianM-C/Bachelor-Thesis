#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate


def histogram(data, bins, fname, label=None, show=False, yscale='', xscale='',
              ylabel='', xlabel='', stacked=False, normed=False, weights=None,
              cumulative=False, use_wigner=False, use_poisson=False):
    fig, ax = plt.subplots()
    ax.hist(data, bins=bins, label=label, stacked=stacked, normed=normed,
            cumulative=cumulative, weights=weights)
    if use_wigner:
        # Wigner distribution
        def wigner_dist(s): return np.pi / 2 * s * np.exp(- np.pi / 4 * s**2)
        bin_size = 1 / 4
        bins = np.arange(0, 4, bin_size)
        # The area of a bar should be the integral of the Wigner distribution
        # between the edges of the bar
        # The are of a bar: A = h * bin_size => h = A / bin_size
        wigner_hist = [1 / bin_size *
                       integrate.quad(wigner_dist, bins[i-1], bins[i])[0]
                       for i in range(1, bins.size)]
        x = np.linspace(0, 4, 100)
        ax.bar(bins[:-1], wigner_hist, width=bin_size,
               align='edge', label='Wigner bar', fill=False, linestyle='-')
        ax.plot(x, wigner_dist(x), 'r-.', label='Wigner')
    if use_poisson:
        def poisson(s): return np.exp(-s)
        x = np.linspace(0, 4, 100)
        ax.plot(x, poisson(x), 'c:', label='Poisson')
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.legend()
    if yscale is 'log':
        ax.set_yscale('log', nonposy='clip')
    if xscale is 'log':
        ax.set_xscale('log')
    fig.savefig(fname, dpi=200)
    if show:
        plt.show()
    plt.close()


def bar_plot(data, fname, label=None, ylabel='', xlabel='', show=False,
             yscale='', xscale='', figsize=None, dpi=None,  axhline_y=0,
             bbox_inches=None):
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(range(data.size), data, label=label)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    if axhline_y:
        ax.axhline(y=axhline_y)
    ax.legend()
    if yscale is 'log':
        ax.set_yscale('log', nonposy='clip')
    if xscale is 'log':
        ax.set_xscale('log')
    fig.savefig(fname, dpi=dpi, bbox_inches=bbox_inches)
    if show:
        plt.show()
    plt.close()
