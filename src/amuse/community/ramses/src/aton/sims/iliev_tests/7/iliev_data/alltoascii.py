#!/usr/bin/python

import os
import glob

datasets = []
for filename in glob.glob('*_add*.bin'):
    i = filename.find('_add')
    j = filename.find('.bin')
    name = filename[:i]
    label = filename[i+4:j]
    datasets.append((name, label))

for name, label in datasets:
    file1 = '%s%s.bin' % (name, label)
    file2 = '%s_add%s.bin' % (name, label)
    out_profile = '%s%s.txt' % (name, label)
    out_slice = '%s_slice%s.txt' % (name, label)
    command = 'python toascii.py "%s" "%s" "%s" "%s"' % (
        file1, file2, out_profile, out_slice)
    print command
    os.system(command)
