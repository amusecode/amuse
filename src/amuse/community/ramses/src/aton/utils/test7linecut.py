#!/usr/bin/python

import glob
import math
import os
import sys

def load_info(output_dir):
    info = {}
    info_path = glob.glob(os.path.join(output_dir, 'info_*'))[0]
    for line in file(info_path):
        parts = line.strip().split()
        if len(parts) == 3 and parts[1] == '=':
            info[parts[0]] = parts[2]
    return info

def calc_units(info):
    mH = 1.66e-24
    kB = 1.38062e-16
    scale_d = float(info['unit_d'])
    scale_t = float(info['unit_t'])
    scale_l = float(info['unit_l'])
    scale_v = scale_l / scale_t
    scale_T2 = mH/kB * scale_v**2
    scale_pressure = scale_d * scale_v**2
    return {'rho': scale_d / mH,
            'p': scale_pressure,
            'x': 1.0,
            't2': scale_T2,
            }

info = load_info('output_00002')
units = calc_units(info)
N = 2**int(info['levelmin'])

points = []

for line in sys.stdin:
    x, y, z, dx, icpu, ilevel, v1, v2, v3, v4, v5, v6 = map(float, line.split())

    mid = (0.5*N + 0.5)/N
    d = math.sqrt((y-mid)**2 + (z-mid)**2)
    if d > 0.5/N:
        continue

    # The vi is related to uold(i) as follows:
    #   v1 = uold(1)
    #   v2 = uold(2) / uold(1)
    #   v3 = uold(3) / uold(1)
    #   v4 = uold(4) / uold(1)
    #   v5 = (gamma-1) * (uold(5) - ekk)
    #   v6 = uold(6) / uold(1)
    # See hydro/output_hydro.f90.

    rho = v1 * units['rho']
    xion = v6 * units['x']
    xneutral = 1.0 - xion
    p = v5 * units['p']
    T2 = v5/v1 * units['t2']
    temperature = T2 / (1 + xion)

    vel = math.sqrt(v2**2 + v3**2 + v4**2)
    c = math.sqrt(1.66667 * v5 / v1)
    mach = vel / c

    points.append((x, rho, xneutral, p, temperature, mach, xion))

points.sort()
for p in points:
    print ' '.join(map(str, p))
