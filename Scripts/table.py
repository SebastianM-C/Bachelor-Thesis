import os
import numpy as np


def cell(x):
    if x == 0:
        return "\\cellcolor{red}" + str(x)
    in_rebde = np.in1d(x, rebde)
    if any(in_rebde):
        return "\\cellcolor{yellow}" + str(x)
    in_reuna = np.in1d(x, reuna)
    if any(in_reuna):
        return "\\cellcolor{blue!50}" + str(x)
    in_reuns = np.in1d(x, reuns)
    if any(in_reuns):
        return "\\cellcolor{green!80}" + str(x)
    return str(x)


A = 1
B = [0, 0.1 * A, 0.4 * A, 0.5 * A]
D = [0, 0.4 * A, 0.5 * A]
N = [40, 60, 80, 100, 120, 140, 160]

# Input parameters
b = 2
d = 2
n = 3

os.chdir("../Output/B" + str(B[b]) + " D" + str(D[d]) + " N" + str(N[n]))

E = np.loadtxt("hamilt.dat", usecols=(1,), unpack=True)
rebde = np.loadtxt("rebde.dat", usecols=(1,), unpack=True)
reuna = np.loadtxt("reuna.dat", usecols=(1,), unpack=True)
reuns = np.loadtxt("reuns.dat", usecols=(1,), unpack=True)

out = "table.tex"

p_rebde = np.pad(rebde, pad_width=(0, E.size - rebde.size), mode='constant')
p_reuna = np.pad(reuna, pad_width=(0, E.size - reuna.size), mode='constant')
p_reuns = np.pad(reuns, pad_width=(0, E.size - reuns.size), mode='constant')

files = [E, p_rebde, p_reuna, p_reuns]

with open(out, "w") as f:
    f.write("\\documentclass{article}\n\n")
    f.write("\\usepackage[margin=0.2in]{geometry}\n\\usepackage{longtable}\n")
    f.write("\\usepackage[table]{xcolor}\n\n")
    f.write("\\begin{document}\n\n")
    f.write("\\begin{longtable}{" + " | ".join(["c"] * 4) + "}\n")
    f.write("hamilt.dat & rebde.dat & reuna.dat & reuns.dat\t\\\\\n")
    f.write("\\hline\n\\endfirsthead\n")
    for row in np.nditer(files, order='C'):
        f.write('\t' + " & ".join(cell(x) for x in row) + " \\\\\n")
    f.write("\\end{longtable}")
    f.write("\n\n\\end{document}")

os.chdir("../../Scripts")
