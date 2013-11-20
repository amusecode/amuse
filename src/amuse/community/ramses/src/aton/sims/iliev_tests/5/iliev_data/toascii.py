#!/usr/bin/python

import struct
import sys
import math

file1 = sys.argv[1]
file2 = sys.argv[2]

s1 = file(file1).read()
s2 = file(file2).read()

def detect_header_width(s):
    int64_header = struct.unpack('q', s[:8])[0]
    int32_header = struct.unpack('i', s[:4])[0]
    assert int64_header == 12 or int32_header == 12
    if int64_header == 12:
        return 'q', 2
    else:
        return 'i', 1

htype, width = detect_header_width(s1)

i = struct.calcsize(htype + 'iii' + htype)
header_size, n1, n2, n3, ignore = struct.unpack(htype + 'iii' + htype, s1[:i])
s1 = s1[i:]
s2 = s2[i:]

assert header_size == 12
assert n1 == 128
assert n2 == 128
assert n3 == 128

float_size = struct.calcsize('f')
array_stride = (2*width + n1*n2*n3) * float_size

print '# x, rho, xneutral, pressure, temperature, mach, xion'

bins = [[0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for i in range(128)]

for i in range(n1):
    for j in range(n2):
        for k in range(n3):
            bin = int(math.sqrt(i*i + j*j + k*k))
            if bin >= len(bins):
                continue
            
            array_index = (width + i + j*n1 + k*n1*n2) * float_size

            q1s, q2s = [], []
            for n in range(3):
                pos = array_index + n*array_stride
                q1 = struct.unpack_from('f', s1, offset=pos)[0]
                q2 = struct.unpack_from('f', s2, offset=pos)[0]
                q1s.append(q1)
                q2s.append(q2)
            
            rho = q2s[0]
            xneutral = q1s[0]
            pressure = q1s[1]
            temperature = q1s[2]
            mach = q2s[1]
            xion = q2s[2]

            bins[bin][0] += 1
            bins[bin][1] += rho
            bins[bin][2] += xneutral
            bins[bin][3] += pressure
            bins[bin][4] += temperature
            bins[bin][5] += mach
            bins[bin][6] += xion

for i in range(len(bins)):
    r = float(i) / 128.0
    v = bins[i][1:]
    for j in range(len(v)):
        v[j] = str(v[j] / bins[i][0])
    print r, ' '.join(v)
