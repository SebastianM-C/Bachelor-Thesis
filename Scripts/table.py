import os
import numpy as np
from timeit import default_timer as timer


def color(x):
    if x == 1:
        return "\\cellcolor{yellow}"
    if x == 2:
        return "\\cellcolor{blue!50}"
    if x == 3:
        return "\\cellcolor{green!80}"
    return ''


A = 1
B = [0, 0.1 * A, 0.4 * A, 0.5 * A]
D = [0, 0.4 * A, 0.5 * A]
N = [40, 60, 80, 100, 120, 140, 160]

# Input parameters
b = 2
d = 2
n = 4

start = timer()
os.chdir("../Output/B" + str(B[b]) + " D" + str(D[d]) + " N" + str(N[n]))

# Read data
E = np.loadtxt("hamilt.dat", usecols=(1,), unpack=True)
rebde = np.loadtxt("rebde.dat", usecols=(1,), unpack=True)
reuna = np.loadtxt("reuna.dat", usecols=(1,), unpack=True)
reuns = np.loadtxt("reuns.dat", usecols=(1,), unpack=True)

# Intersections
in_rebde = np.in1d(E, rebde, assume_unique=True)
in_reuna = np.in1d(E, reuna, assume_unique=True)
in_reuns = np.in1d(E, reuns, assume_unique=True)

# Padding
p_rebde = np.pad(rebde, pad_width=(0, E.size - rebde.size), mode='constant')
p_reuna = np.pad(reuna, pad_width=(0, E.size - reuna.size), mode='constant')
p_reuns = np.pad(reuns, pad_width=(0, E.size - reuns.size), mode='constant')

files = np.array([E, p_rebde, p_reuna, p_reuns, in_rebde, in_reuna, in_reuns])

with open("table.tex", "w") as f:
    f.write("\\documentclass{article}\n\n")
    f.write("\\usepackage[margin=0.2in]{geometry}\n\\usepackage{longtable}\n")
    f.write("\\usepackage[table]{xcolor}\n\n")
    f.write("\\begin{document}\n\n")
    f.write("\\begin{longtable}{" + " | ".join(["c"] * 4) + "}\n")
    f.write("hamilt.dat & rebde.dat & reuna.dat & reuns.dat\t\\\\\n")
    f.write("\\hline\n\\endfirsthead\n")
    for row in range(files.shape[1]):
        # Find the representation of the energy level
        c = 0
        for x in range(1, 4):
            if files[x + 3][row]:
                c = x

        line = color(c) + str(files[0][row]) + " & " + \
            " & ".join(color(x) + str(files[x][row]) for x in range(1, 4))
        f.write('\t' + line + " \\\\\n")
    f.write("\\end{longtable}")
    f.write("\n\n\\end{document}")

os.chdir("../../Scripts")
end = timer()
print(end - start)
