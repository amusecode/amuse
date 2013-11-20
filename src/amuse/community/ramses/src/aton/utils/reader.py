#!/usr/bin/python

import glob
import os
import os.path
import sys

prefix = ''
if sys.argv[0]:
    prefix = os.path.dirname(sys.argv[0])
amr2map_bin = os.path.join(prefix, '../f90/amr2map')
part2map_bin = os.path.join(prefix, '../f90/part2map')
histo_bin = os.path.join(prefix, '../f90/histo')
getstarlist_bin = os.path.join(prefix, '../f90/getstarlist')

def load_info(output_dir):
    info = {}
    files = glob.glob(os.path.join(output_dir, 'info_*'))
    if not files:
        return None
    info_path = files[0]
    for line in file(info_path):
        parts = line.strip().split()
        if len(parts) == 3 and parts[1] == '=':
            info[parts[0]] = parts[2]
    info['aexp'] = float(info['aexp'])
    info['z'] = 1.0/info['aexp'] - 1.0

    unit_l = info['unit_l'] = float(info['unit_l'])
    unit_d = info['unit_d'] = float(info['unit_d'])
    unit_t = info['unit_t'] = float(info['unit_t'])
    unit_v = unit_l / unit_t
    mH = 1.6600000e-24
    kB = 1.3806200e-16
    info['unit_T2'] = mH/kB * unit_v**2

    return info

def load_radiation(output_dir):
    def parse_shard(f):
        d = {}
        section = d.get('default', {})
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('#'):
                name = line[1:].strip().lower()
                section = d.setdefault(name, {})
            else:
                i = line.find(' ')
                label = line[:i]
                value = line[i+1:]
                section[label] = value
        return d

    avg = {}
    n = {}
    rad_shards = []
    for path in glob.glob(os.path.join(output_dir, 'radstats_*')):
        rad_shard = parse_shard(file(path))
        rad_shards.append(rad_shard)
        s = rad_shard['statistics']
        for k in s:
            v = float(s[k].split()[0])
            avg[k] = avg.get(k, 0.0) + v
            n[k] = n.get(k, 0) + 1
    for k in avg:
        avg[k] = avg[k] / n[k]
    return avg, rad_shards

def exists(path):
    try:
        s = os.stat(path)
        return True
    except OSError, e:
        return False

def amr2map(output_dir, out_txt, type):
    if exists(out_txt):
        return
    amr2map_params = {
        'inp': output_dir,
        'out': out_txt,
        'fil': 'ascii',
        'typ': str(type)  # ionized fraction
        }
    amr2map_command = amr2map_bin
    for param in amr2map_params.items():
        amr2map_command += ' -%s %s' % param
    print >>sys.stderr, 'running:', amr2map_command
    assert os.system(amr2map_command) == 0

def part2map(output_dir, out_txt, stars):
    if exists(out_txt):
        return
    part2map_params = {
        'inp': output_dir,
        'out': out_txt,
        'fil': 'ascii',
        'str': stars and 'true' or 'false',
        }
    part2map_command = part2map_bin
    for param in part2map_params.items():
        part2map_command += ' -%s %s' % param
    print >>sys.stderr, 'running:', part2map_command
    assert os.system(part2map_command) == 0

def histo(output_dir, out_txt, extra_params={}):
    if exists(out_txt):
        return
    histo_params = {
        'inp': output_dir,
        'out': out_txt,
        'fil': 'ascii',
        }
    histo_params.update(extra_params)
    histo_command = histo_bin
    for param in histo_params.items():
        histo_command += ' -%s %s' % param
        print >>sys.stderr, 'running:', histo_command
    assert os.system(histo_command) == 0

def getstarlist(output_dir, out_txt):
    if exists(out_txt):
        return
    getstarlist_params = {
        'inp': output_dir,
        'out': out_txt,
        }
    getstarlist_command = getstarlist_bin
    for param in getstarlist_params.items():
        getstarlist_command += ' -%s %s' % param
    print >>sys.stderr, 'running:', getstarlist_command
    assert os.system(getstarlist_command) == 0

def gnuplot(script):
    tmp_gp = 'tmp.gp'
    file(tmp_gp, 'w').write(script)
    gnuplot_command = 'gnuplot %s' % tmp_gp
    print >>sys.stderr, 'running:', gnuplot_command
    return os.system(gnuplot_command)

def scan_outputs(prefix):
    outputs = glob.glob(prefix + '/output_*')
    outputs.sort()
    for output_dir in outputs:
        dir = '/output_'
        j = output_dir.rfind('/output_')
        if j == -1:
            print 'skipping', output_dir
            continue
        i = int(output_dir[j + len(dir):])
        info = load_info(output_dir)
        if not info:
            print 'skipping', output_dir
            continue
        yield i, output_dir, info
