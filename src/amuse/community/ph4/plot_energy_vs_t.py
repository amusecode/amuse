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
    
    #if infile == '' or quantity == '':
        #print 'plotq -f file -q quantity -t linetype'
        #sys.exit(0)
    
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
                elif cols[1] == 'Multiple':
                    bin_id = int(cols[2])
                    #if log: q = log10(q)
                    strsplit = cols[10].split('=');
                    if len(strsplit) == 2:
                        mykT = float(strsplit[1])
                        mydata = [tlist[len(tlist)-1], bin_id, mykT]
                        qlist.append(mydata)
                #elif re.search(quantity+'\[.*\]'+'=', cols[1]):
                    #s = re.search(quantity+'\[.*\]'+'=', cols[1]).group()
                    # Messy! There must surely be a better way to do this...
                    #i1 = string.index(s,'[')
                    #i2 = string.index(s,']')
                    #nq = int(s[i1+1:i2])
                    #clist = []
                    #for i in range(2,2+nq):
                        #q = float(cols[i])
                        #if log: q = log10(q)
                        #clist.append(q)
                    #qlist.append(clist)
    f.close()
    
    multiple_index = []
    i=0
    while i<len(qlist):
      j=0
      original = True
      while j<len(multiple_index):
        if multiple_index[j] == qlist[i][1]:
          original = False 
        j += 1
      if original:
        multiple_index.append(qlist[i][1])
      i += 1
    
    counter = 0
    i=0
    while i<len(multiple_index):
      myidata = []
      myxdata = []
      myydata = []
      j=0
      while j<len(qlist):
        if multiple_index[i] == qlist[j][1]:
          myidata.append(multiple_index[i])
          myxdata.append(qlist[j][0])
          myydata.append(qlist[j][2])
          if qlist[j][0] == 300.0:
            counter += 1
        j += 1
      plot(myxdata, myydata)
      hold('on')
      i += 1
    
    show()
    
    
    i=0
    while i<len(myxdata):
      print myxdata[i]
    
      i += 1
    print counter
    
    """
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
    """
