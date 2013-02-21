#!/usr/bin/env python

import getopt
import string
from pylab import *
import re
import numpy
from math import log10

# Display a time sequence from a set of (not yet standard) AMUSE log
# files.  Files, colors, and styles should be separated by ',' with no
# spaces.

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:f:lq:o:t:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)
    
    infiles = ['']
    quantity = ''
    log = 0
    outfile = ''

    # Separate the point colors and the line styles.  Replicate the
    # last color/style in case of multiple files.
    
    colors = ['b','g','r','c','m','k']
    ls = ['-']

    for o, a in opts:
        if o == "-c":
            colors = a.split(',')
        elif o == "-f":
            infiles = a.split(',')
        elif o == "-l":
            log = 1
        elif o == "-o":
            outfile = a
        elif o == "-q":
            quantity = a
        elif o == "-t":
            ls = a.split(',')
        else:
            print "unexpected argument", o
    
    if infiles[0] == '' or quantity == '':
        print 'plotq -f file-list -q quantity -t linetype'
        sys.exit(0)

    li = len(infiles)
    lc = len(colors)
    if lc < li:
        for i in range(li-lc):
            colors.append(colors[lc-1])
    lt = len(ls)
    if lt < li:
        for i in range(li-lt):
            ls.append(ls[lt-1])

    #print infiles
    #print colors
    #print ls

    i = 0
    for infile in infiles:

        tlist = []
        qlist = []
        nq = 1
    
        f = open(infile, 'r')
        for line in f:
            if len(line) > 0:
                cols = line.split()
                if len(cols) > 0 and cols[0] == '%%%':
                    if cols[1] == 'time=':
                        tlist.append(float(cols[2]))
                    elif cols[1] == quantity+'=':
                        q = float(cols[2])
                        if log: q = log10(q)
                        qlist.append(q)
                    elif re.search(quantity+'\[.*\]'+'=', cols[1]):
                        s = re.search(quantity+'\[.*\]'+'=', cols[1]).group()
                        # Messy! There must surely be a better way to do this...
                        i1 = string.index(s,'[')
                        i2 = string.index(s,']')
                        nq = int(s[i1+1:i2])
                        clist = []
                        for i in range(2,2+nq):
                            q = float(cols[i])
                            if log: q = log10(q)
                            clist.append(q)
                        qlist.append(clist)
        f.close()
        
        if len(qlist) > 0:
            if log: quantity = 'log10 '+quantity
            if nq == 1:
                plot(tlist, qlist, colors[i]+ls[i])
            else:
                qplot = numpy.array(qlist)
                for i in range(nq):
                    plot(tlist, qplot[:,i], lt, linewidth=1.0)
            xlabel('time')
            ylabel(quantity)
        else:
            print 'No data to plot in '+infile+'.'

        i += 1

    # Save a file if specified; otherwise, display to the screen.

    if not outfile == '':
        if outfile == '.': outfile = quantity+'.pdf'
        outfile = outfile.replace('/','_')
        print 'saving to', outfile
        savefig(outfile)
    else:
        show()
