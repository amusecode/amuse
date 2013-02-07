#!/usr/bin/env python

import getopt
import string
from pylab import *
import re
import numpy
from math import log10

# Display a time sequence from a (not yet standard) AMUSE log file.

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:lq:t:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)
    
    infile = ''
    quantity = ''
    log = 0
    lt = 'b-'
    
    for o, a in opts:
        if o == "-f":
            infile = a
        elif o == "-l":
            log = 1
        elif o == "-q":
            quantity = a
        elif o == "-t":
            lt = a
        else:
            print "unexpected argument", o
    
    if infile == '' or quantity == '':
        print 'plotq -f file -q quantity -t linetype'
        sys.exit(0)
    
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
            plot(tlist, qlist, lt)
        else:
            qplot = numpy.array(qlist)
            for i in range(nq):
                plot(tlist, qplot[:,i], lt, linewidth=1.0)
        xlabel('time')
        ylabel(quantity)
        show()
    else:
        print 'No data to plot.'
