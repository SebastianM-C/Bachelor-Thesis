#!/usr/bin/env python

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from custom_parser import parse

from tools import find, get_input


def main(N):
    files = find('alpha.txt', '../Output')
    values = []
    avg_deltaE = 0
    count = 0
    for f in files:
        b, d, n = get_input(f)
        if int(n) == N:
            values.append([b, np.asscalar(np.loadtxt(f))])
            avg_deltaE += np.loadtxt(f.replace('alpha', 'stable'))[1]
            count += 1

    avg_deltaE = avg_deltaE / count
    print('Found', count, 'values\navg deltaE', avg_deltaE)

    # Plot the results
    fig, ax = plt.subplots()
    ax.scatter(np.array(values)[:, 0], np.array(values)[:, 1],
               label='$\\Delta E\\approx' + '{:.3f}'.format(avg_deltaE) + '$')
    ax.set_ylabel('$\\alpha$')
    ax.set_xlabel('B')
    ax.legend(markerscale=0)
    plt.show()
    fig.savefig('../Statistics/alpha_N' + str(N) + '_dE' +
                '{:.3f}'.format(avg_deltaE) + '.png', dpi=200)
    plt.close()


if __name__ == '__main__':
    B, D, N = parse()

    for n in N:
        main(n)
