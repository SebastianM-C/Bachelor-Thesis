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
              count=16, fit=False, max_e=0):
    fig, ax = plt.subplots()
    n, bin_edges, _ = ax.hist(data, bins=bins, label=label, stacked=stacked,
                              normed=normed, cumulative=cumulative,
                              weights=weights)

    bin_size = 4 / (count - 1)
    bins = np.linspace(0, 4, count)
    x = np.linspace(0, 4, 100)

    if use_wigner:
        if cumulative:
            wigner_hist = [integrate.quad(wigner, 0, bins[i])[0]
                           for i in range(1, bins.size)]

            def wigner_dist_int(s): return 1 - np.exp(- np.pi / 4 * s**2)
            ax.plot(x, wigner_dist_int(x), 'r-.',
                    label=r'$\int_0^s P_{W}(x)\mathrm{d}x$')
            ax.bar(bins[:-1], wigner_hist, width=bin_size,
                   align='edge', fill=False, linestyle='-',
                   label=r'$\sum_0^S P_W(x)$')
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
        ax.plot(x, poisson(x), 'c:', label='Poisson')
        poisson_hist = [1 / bin_size *
                        integrate.quad(poisson, bins[i-1], bins[i])[0]
                        for i in range(1, bins.size)]
        ax.bar(bins[:-1], poisson_hist, width=bin_size, align='edge',
               fill=False, linestyle=(':'), edgecolor='magenta',
               label='Poisson bar')

    if fit:
        # Fit the data with a superposition between the Poisson and Wigner
        # distributions
        x_data = np.linspace(0, 4, count - 1)
        y_data = n[0] + n[1] + n[2]
        alpha, pcov = curve_fit(model, x_data, y_data, bounds=(0, 1))
        fname_a = 'alpha' + \
            ('_max_e_' + str(max_e) + '.txt' if max_e else '.txt')
        alpha.tofile(fname_a, sep=' ')     # save the value
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
             bbox_inches=None, title=''):
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(range(data.size), data, label=label)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    if axhline_y:
        ax.axhline(y=axhline_y)
    ax.legend(title=title)
    if yscale is 'log':
        ax.set_yscale('log', nonposy='clip')
    if xscale is 'log':
        ax.set_xscale('log')
    fig.savefig(fname, dpi=dpi, bbox_inches=bbox_inches)
    if show:
        plt.show()
    plt.close()
