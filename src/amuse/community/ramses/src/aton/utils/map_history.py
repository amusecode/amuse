import os.path
import reader

def average(output_dir, label):
    try:
        f = file(os.path.join(output_dir, label + '_map.txt'))
    except IOError:
        return 0.0

    s = 0.0
    n = 0
    for line in f:
        line = line.strip()
        if not line: continue
        parts = line.split()
        s += float(parts[2])
        n += 1
    return s / n

f = file('map_history.txt', 'w')
for i, output_dir, info in reader.scan_outputs('.'):
    if i == 1: continue
    print output_dir

    reader.amr2map(output_dir, os.path.join(output_dir, 'temperature_map.txt'), 8)

    avg_xion_volume = average(output_dir, 'xion')
    if avg_xion_volume is None:
        continue
    avg_xion_mass = 0.0
    avg_density = average(output_dir, 'density') * info['unit_d']
    avg_T2 = average(output_dir, 'temperature') * info['unit_T2']

    z = info['z']

    print >>f, z, avg_xion_volume, avg_xion_mass, avg_T2, avg_density
    f.flush()
