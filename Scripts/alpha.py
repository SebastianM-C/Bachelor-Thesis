#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from custom_parser import parse

from tools import find, get_input


def main(N):
    files = find('alpha.txt', '../Output')
    fig, ax = plt.subplots()
    msize = 8
    for n_i in N:
        values = []
        count = 0
        avg_deltaE = 0
        for f in files:
            b, d, n = get_input(f)
            if int(n) == n_i:
                values.append([b, np.asscalar(np.loadtxt(f))])
                avg_deltaE += np.loadtxt(f.replace('alpha', 'stable'))[1]
                count += 1
        avg_deltaE = avg_deltaE / count
        print('Found', count, 'values for N =', n_i,
              '\navg deltaE', avg_deltaE)
        # Plot the results
        ax.plot(np.array(values)[:, 0], np.array(values)[:, 1], 'o',
                label=r'$\Delta E\approx' + '{:.3f}'.format(avg_deltaE) +
                '$', markersize=msize)
        msize -= 2
    ax.set_ylabel('$\\alpha$')
    ax.set_xlabel('$B$')
    ax.set_ylim([0, 1.1])
    ax.set_xlim([0, 1])
    ax.legend()
    # plt.show()
    fig.savefig('../Statistics/alpha_N' + str(N) + '_dE' +
                '{:.3f}'.format(avg_deltaE) + '.png', dpi=200)
    plt.close()


if __name__ == '__main__':
    N = parse(n_only=True)
    main(N)
