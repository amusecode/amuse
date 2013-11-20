#!/usr/bin/python

"""
Reads the output from amr2cell and outputs the averaged quantities in spherical
shells about the centre of the cube.
"""

import math
import sys

N = 256
bins = [[0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for i in range(N)]

middle = 0.5

# Units
# These need to be kept in sync with amr/units.f90!
mH = 1.66e-24
kB = 1.38062e-16
scale_d = mH * 1.0e-3
scale_t = 3.156e15
scale_l = 2*4.629e22
scale_v = scale_l / scale_t
scale_T2 = mH/kB * scale_v**2
scale_pressure = scale_d * scale_v**2

for line in sys.stdin:
    x, y, z, dx, icpu, ilevel, v1, v2, v3, v4, v5, v6 = map(float, line.split())

    # The vi is related to uold(i) as follows:
    #   v1 = uold(1)
    #   v2 = uold(2) / uold(1)
    #   v3 = uold(3) / uold(1)
    #   v4 = uold(4) / uold(1)
    #   v5 = (gamma-1) * (uold(5) - ekk)
    #   v6 = uold(6) / uold(1)
    # See hydro/output_hydro.f90.

    rho = v1 * scale_d / mH
    xion = v6
    xneutral = 1.0 - v6
    p = v5 * scale_pressure
    T2 = v5/v1 * scale_T2
    temperature = T2 / (1 + xion)

    vel = math.sqrt(v2**2 + v3**2 + v4**2)
    c = math.sqrt(1.66667 * v5 / v1)
    mach = vel / c

    v = [1, rho, xneutral, p, temperature, mach, xion]
    r = 2*math.sqrt((x-middle)**2 + (y-middle)**2 + (z-middle)**2)
    if r >= 1.0:
        continue
    i = int(r*N)

    for j in range(len(v)):
        bins[i][j] += v[j]

for i in range(N):
    if bins[i][0] == 0:
        continue

    for j in range(1,len(bins[i])):
        bins[i][j] /= bins[i][0]

    r = float(i)/N
    print r, ' '.join(map(str, bins[i][1:]))
