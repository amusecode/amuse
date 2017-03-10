#!/usr/bin/env python

import sys
from math import log10, sqrt
import numpy as np
import matplotlib.pyplot as plt
import string
import re

# Display a time sequence from a set of (not yet standard) AMUSE log
# files.  Files, colors, and styles should be separated by spaces.

if __name__ == '__main__':

    # Command-line arguments and defaults:

    infiles = []				# -f [list]
    log = 0					# -L (sets log = 1)
    outfile = ''				# -o filename or .
    quantities = []				# -q [list]
    layout = ''					# -l (rows)x(cols)
    space = 0.3					# -S space

    # Separate the point colors and the line styles.  Replicate the
    # last color/style in case of multiple files.

    colors = ['b','g','r','c','m','k']		# -c [list]
    ls = ['-']					# -s [list]

    # Parse the argument list the old-fashioned way to enable multiple
    # arguments (unaccountably not supported by getopt):

    i = 1
    n = len(sys.argv)
    while i < n:
        if sys.argv[i][0] == '-':

            if len(sys.argv[i]) > 1:

                if sys.argv[i][1] == 'c':

                    colors = []

                    # Accumulate colors until we find a leading '-' or
                    # reach the end of the list.

                    while i < n-1 and sys.argv[i+1][0] != '-':
                        i += 1
                        colors.append(sys.argv[i])

                elif sys.argv[i][1] == 'f':

                    # Accumulate names until we find a leading '-' or
                    # reach the end of the list.

                    while i < n-1 and sys.argv[i+1][0] != '-':
                        i += 1
                        infiles.append(sys.argv[i])

                elif sys.argv[i][1] == 'L':

                    log = 1

                elif sys.argv[i][1] == 'l':

                    if i < n-1 and sys.argv[i+1][0] != '-':
                        i += 1
                        layout = sys.argv[i]

                elif sys.argv[i][1] == 'o':

                    if i < n-1 and sys.argv[i+1][0] != '-':
                        i += 1
                        outfile = sys.argv[i]

                elif sys.argv[i][1] == 'q':

                    # Accumulate quantities until we find a leading
                    # '-' or reach the end of the list.

                    while i < n-1 and sys.argv[i+1][0] != '-':
                        i += 1
                        quantities.append(sys.argv[i])

                elif sys.argv[i][1] == 'S':

                    if i < n-1 and sys.argv[i+1][0] != '-':
                        i += 1
                        space = float(sys.argv[i])

                elif sys.argv[i][1] == 's':

                    ls = []

                    # Accumulate line styles until we find a leading
                    # '-' or reach the end of the list.  Problem: '-'
                    # is a legal line style, so accept a lone '-' in
                    # this case.

                    while i < n-1 \
                            and (sys.argv[i+1][0] != '-'
                                  or len(sys.argv[i+1]) == 1):
                        i += 1
                        ls.append(sys.argv[i])

                else:
                    print 'unknown option', sys.argv[i]

        i += 1
    
    li = len(infiles)
    lq = len(quantities)
    if li == 0 or lq == 0:
        print 'plotq -f file-list -q quantity-list'
        sys.exit(0)

    # Extend the color list if neessary.

    lc = len(colors)
    if lc < li:
        for i in range(li-lc):
            colors.append(colors[lc-1])

    # Extend the line style list if necessary.

    lt = len(ls)
    if lt < li:
        for i in range(li-lt):
            ls.append(ls[lt-1])

    plt.figure(figsize=(12,8))
    plt.subplots_adjust(wspace=space, hspace=space)

    # Determine the plot layout (mxn).

    if layout == '':

        # Automatic.

        m = int(sqrt(float(lq)))
        if m*m < lq: m += 1
        n = lq/m
        if m*n < lq: n += 1

    else:

        # User-specified:

        m,n = layout.split('x')
        m = int(m)
        n = int(n)

    i = 100*m+10*n				# assume m, n < 10
    figs = []
    for ii in range(lq):
        figs.append(plt.subplot(i+ii+1))

    print 'infiles:    ',
    for f in infiles: print f,
    print ''
    print 'quantities: ',
    for q in quantities: print q,
    print ''
    print 'colors:     ',
    for ic in range(li): print colors[ic],
    print ''
    print 'styles:     ',
    for i in range(li): print ls[i],
    print ''
    print m, 'x', n, 'plot layout'

    ifile = 0
    for infile in infiles:			# loop over input files

        # Extract data from infile.

        tlist = []
        qlist = [[] for iq in range(lq)]	# list of lq empty lists
        nq = [1 for iq in range(lq)]		# default qlist length

        f = open(infile, 'r')
        for line in f:				# process data line by line
            if len(line) > 0:
                cols = line.split()
                if len(cols) > 0 and cols[0] == '%%%':

                    if cols[1] == 'time=':
                        tlist.append(float(cols[2]))

                    else:

                        foundq = False
                        for iq in range(lq):
                            if cols[1] == quantities[iq]+'=':
                                foundq = True
                                q = float(cols[2])
                                if log: q = log10(q)
                                qlist[iq].append(q)

                        if not foundq:
                            for iq in range(lq):
                                qq = quantities[iq]
                                if re.search(qq+'\[.*\]'+'=', cols[1]):
                                    s = re.search(qq+'\[.*\]'+'=',
                                                   cols[1]).group()
                                    # Messy! There must surely be a
                                    # better way to do this...
                                    i1 = string.index(s,'[')
                                    i2 = string.index(s,']')
                                    nq[iq] = int(s[i1+1:i2])
                                    clist = []
                                    for i in range(2,2+nq[iq]):
                                        q = float(cols[i])
                                        if log: q = log10(q)
                                        clist.append(q)
                                    qlist[iq].append(clist)

        f.close()

        # Plot the data.

        for iq in range(lq):
            qq = quantities[iq]
            if len(qlist[iq]) > 0:
                if log: qq = 'log10 '+qq
                if nq[iq] == 1:
                    figs[iq].plot(tlist, qlist[iq], colors[ifile]+ls[ifile])
                else:
                    qplot = np.array(qlist[iq])
                    for i in range(nq[iq]):
                        figs[iq].plot(tlist, qplot[:,i],
                             colors[ifile]+ls[ifile], linewidth=1.0)
                figs[iq].set_xlabel('time')
                figs[iq].set_ylabel(qq)
            else:
                print 'No', qq, 'data to plot in '+infile+'.'

        ifile += 1

    # Save to a file if specified; otherwise, display to the screen.

    if outfile == '':
        plt.show()
    else:
        if outfile == '.':
            qq = ''
            for iq in range(lq):
                if iq > 0: qq += '+'
                qq += quantities[iq]
            outfile = qq+'.pdf'
        outfile = outfile.replace('/','_')
        print 'saving to', outfile
        plt.savefig(outfile)
