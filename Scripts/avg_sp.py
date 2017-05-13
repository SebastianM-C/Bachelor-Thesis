#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from custom_parser import parse

from tools import find, get_input


def main(B, D, N):
    files = find('avg_sp.txt', '../Output')
    fig, ax = plt.subplots()
    msize = 8
    for n_i in N:
        values = []
        count = 0
        for f in files:
            b, d, n = get_input(f)
            if int(n) == n_i:
                values.append([b, np.loadtxt(f)])
                count += 1
        print('Found', count, 'values for N =', n_i)
        results = np.zeros(count, [('b', np.float64),
                                   ('avg_sp', np.float64, 3)])

        for i in range(count):
            results[i] = tuple(values[i])
        # Plot the results
        ax.plot(results['b'], results['avg_sp'][:, 0], '^',
                label='reuna, $N=' + str(n_i) + '$', markersize=msize)
        ax.plot(results['b'], results['avg_sp'][:, 1], 'v',
                label='reuns, $N=' + str(n_i) + '$', markersize=msize)
        ax.plot(results['b'], results['avg_sp'][:, 2], 'o',
                label='rebde, $N=' + str(n_i) + '$', markersize=msize)
        msize -= 2
    ax.set_ylabel('$\\frac{E_n^r - E_0^r}{N^r}$')
    ax.set_xlabel('$B$')
    ax.set_xlim([0, 1])
    ax.legend()
    # plt.show()
    fig.savefig('../Statistics/avg_sp_N' + str(N) + '.png', dpi=200)
    plt.close()


if __name__ == '__main__':
    B, D, N = parse()
    main(B, D, N)
