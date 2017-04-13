#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np


def histogram(data, bins, fname, label=None, show=False, yscale='', xscale='',
              ylabel='', xlabel='', stacked=False, normed=False,
              cumulative=False, weights=None):
    fig, ax = plt.subplots()
    ax.hist(data, bins=bins, label=label, stacked=stacked, normed=normed,
            cumulative=cumulative, weights=weights)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.legend()
    if yscale is 'log':
        ax.set_yscale('log', nonposy='clip')
    if xscale is 'log':
        ax.set_xscale('log')
    fig.savefig(fname)
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
