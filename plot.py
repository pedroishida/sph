#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True

xs = np.array([])
ys = np.array([])
rs = np.array([])
ds = np.array([])

with open("sphstar.csv", 'r') as f:
    l, k, m, N, n, dt, h, damp = [float(x) for x in f.readline().split(',')]
    for line in f:
        x, y, d = [float(x) for x in line.split(',')]
        xs = np.append(xs, x)
        ys = np.append(ys, y)
        rs = np.append(rs, (x**2 + y**2)**0.5)
        ds = np.append(ds, d)

m = N * m

dmax = 0.0
for d in ds:
    if d > dmax:
        dmax = d
ds = ds / dmax

radius = (2.0 * k / l)**(1.0 / (2 * (1 + 1.0 / n))) * (m / np.pi)**(1.0 / (2 * (n + 1))) * (n + 1)**0.5
print radius
ro = np.linspace(0, radius, 20)
do = (radius**2 - ro**2)**n

dmax = 0.0
for d in do:
    if d > dmax:
        dmax = d
do = do / dmax

label = "dt = " + str(dt) + "\nh = " + str(h) + "\ndamp = " + str(damp) + "\nN = " + str(N)

plt.plot(rs, ds, 'bo')
plt.plot(ro, do, 'r-')
plt.text(0.1, 0.5, label)
plt.xlabel("r")
plt.ylabel(r"$\frac{\rho}{\rho_{max}}$", fontsize=25, rotation=0, horizontalalignment='right')
plt.ylim([0, 1.1])
plt.show()

plt.scatter(xs, ys, c=ds, marker='o', edgecolors='none')
plt.text(-0.9, -0.9, label)
plt.colorbar(label=r"$\frac{\rho}{\rho_{max}}$")
plt.xlabel("x")
plt.ylabel("y")
plt.show()
