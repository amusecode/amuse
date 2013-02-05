from pylab import *
import sys
import re

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "usage: python binary_hist.py <input filename> <time of histogram>"
        sys.exit(1)
    else:
        fname = sys.argv[1]
        time = float(sys.argv[2])
        f = open(fname, "r")
        inblock = False
        EkTs = []
        for line in f:
            if re.search("%%% time= (\d+\.\d*)", line):
                if float(re.search("%%% time= (\d+\.\d*)", line).group(1)) == time:
                    inblock = True
            if inblock and re.search("%%% Multiple.*E/kT=(-?\d+\.\d+)", line):
                EkTs.append(float(re.search("%%% Multiple.*E/kT=(-?\d+\.\d+)",
                                            line).group(1)))
    
            if inblock and re.search("%%% Emul/E", line):
                inblock = False
        f.close()
    
        if len(EkTs) > 0:
            hist(EkTs)
            show()
        else:
            print "No binaries found at time = %f." % time
    
