#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit


def wigner(s):
    """Wigner distribution"""
    return np.pi / 2 * s * np.exp(- np.pi / 4 * s**2)


def poisson(s):
    """Poisson distribution"""
    return np.exp(-s)


def model(s, alpha):
    """Superposition between Poisson and Wigner distributions"""
    return alpha * poisson(s) + (1 - alpha) * wigner(s)


def histogram(data, bins, fname, label=None, show=False, yscale='', xscale='',
              ylabel='', xlabel='', stacked=False, normed=False, weights=None,
              cumulative=False, use_wigner=False, use_poisson=False, title='',
              bin_size=1/4, fit=False):
    fig, ax = plt.subplots()
    n, bin_edges, _ = ax.hist(data, bins=bins, label=label, stacked=stacked,
                              normed=normed, cumulative=cumulative,
                              weights=weights)
    if use_wigner:
        bins = np.arange(0, 4, bin_size)
        x = np.linspace(0, 4 - bin_size, 100)
        if cumulative:
            wigner_hist = [integrate.quad(wigner, 0, bins[i])[0]
                           for i in range(1, bins.size)]

            def wigner_dist_int(s): return 1 - np.exp(- np.pi / 4 * s**2)
            ax.plot(x, wigner_dist_int(x), 'r-.', label='Wigner')
        else:
            # The area of a bar should be the integral of the Wigner
            # distribution between the edges of the bar
            # The are of a bar: A = h * bin_size => h = A / bin_size
            wigner_hist = [1 / bin_size *
                           integrate.quad(wigner, bins[i-1], bins[i])[0]
                           for i in range(1, bins.size)]
            ax.plot(x, wigner(x), 'r-.', label='Wigner')

        ax.bar(bins[:-1], wigner_hist, width=bin_size,
               align='edge', label='Wigner bar', fill=False, linestyle='-')
    if use_poisson:
        x = np.linspace(0, 4 - bin_size, 100)
        ax.plot(x, poisson(x), 'c:', label='Poisson')

    if fit:
        # Fit the data with a superposition between the Poisson and Wigner
        # distributions
        x = np.linspace(0, 4 - bin_size, 100)
        x_data = np.arange(0, 4 - bin_size, bin_size)
        y_data = n[0] + n[1] + n[2]
        alpha, pcov = curve_fit(model, x_data, y_data, bounds=(0, 1))
        alpha.tofile('alpha.txt', sep=' ')     # save the value
        model_hist = [1 / bin_size *
                      integrate.quad(model, bins[i-1], bins[i], args=alpha)[0]
                      for i in range(1, bins.size)]
        ax.bar(bins[:-1], model_hist, width=bin_size, align='edge', fill=False,
               linestyle='-')
        ax.plot(x, model(x, alpha), 'r',
                label='$\\alpha = ' + '{:.3}'.format(alpha[0]) + '$')
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.legend(title=title)
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
