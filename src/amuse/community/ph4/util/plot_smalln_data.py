import sys
import getopt
import numpy
from scipy.interpolate import interpolate
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animate

# TODO:
#	Automatic determination of plot limits.
#	Automatic determination of dt.
#	Command-line selection of projection axis.
#	Colormap should respond to changes in index array.

def read_data(filename):

    # Read in the data file.  Format is defined in smalln.

    print 'reading data from', filename
    sys.stdout.flush()

    f = open(filename)
    ll = f.readlines()
    f.close()

    print len(ll), 'records read'
    sys.stdout.flush()

    t = []
    id = []
    x = []
    y = []
    z = []
    line = 0

    for l in ll:
        line += 1

	# Format is: time  id1 x1 y1 z1  id2 x2 y2 z2  id3 x3 y3 z3 ...
    
        r = l.split()
        n = (len(r)-1)/4
        if line == 1: print 'n =', n
        if 4*n + 1 != len(r):
            print 'record length error'
            sys.exit(1)

        # Unpack the line.

        tt = float(r[0])
        ii = []
        xx = []
        yy = []
        zz = []
        for i in range(n):
            ii.append(int(r[4*i+1]))
            xx.append(float(r[4*i+2]))
            yy.append(float(r[4*i+3]))
            zz.append(float(r[4*i+4]))

        # Add the data to the master arrays.

        t.append(tt)
        id.append(ii)
        x.append(xx)
        y.append(yy)
        z.append(zz)

    ta = numpy.array(t)
    ia = numpy.array(id)
    xa = numpy.array(x)
    ya = numpy.array(y)
    za = numpy.array(z)

    return len(ta), ta, ia, xa, ya, za

def interpolate_data(t, x, y, dt):
    nt, np = x.shape
    print 'interpolating to dt =', dt
    print 'mean time interval  =', (t[nt-1]-t[0])/nt
    sys.stdout.flush()

    tout = numpy.arange(t[0], t[nt-1], dt)
    nout = len(tout)
    xout = numpy.zeros((nout,np))
    yout = numpy.zeros((nout,np))
    for k in range(np):
        f = interpolate.interp1d(t, x[:,k], kind='cubic', assume_sorted=True)
        xout[:,k] = f(tout)
        f = interpolate.interp1d(t, y[:,k], kind='cubic', assume_sorted=True)
        yout[:,k] = f(tout)

    print 'nout =', nout
    return nout, tout, xout, yout

def animate_data(t, x, y, id, lx, ly, xmax):

    print 'animating data'
    colormap = ['r', 'y', 'b', 'm', 'g', 'c', 'k', 'k', 'k', 'k', 'k', 'k']

    def update(frame):

        # What idiot decided that the syntax for set_offsets should be
        # different from the syntax for scatter??

        off = []
        for j in range(np):
            off.append(x[frame,j])
            off.append(y[frame,j])
        scat.set_offsets(off)
        fig.suptitle('time ='+'%6.2f'%(t[frame]))

        return scat,

    nt,np = x.shape
    fig = plt.figure()
    scat = plt.scatter(x[0,:], y[0,:], c=colormap[:np], s=30)
    #plt.axis('square')	# fails in ubuntu
    plt.xlabel(lx)
    plt.ylabel(ly)
    plt.xlim(-xmax, xmax)
    plt.ylim(-xmax, xmax)
    plt.grid()

    ani = animate.FuncAnimation(fig, update, frames=range(nt), interval=10)
    plt.show()

def plot_data(t, x, lx, ly):
    nt,np = x.shape
    fig = plt.figure()
    plt.plot(t, x, 'b-')
    plt.xlabel(lx)
    plt.ylabel(ly)
    plt.grid()
    plt.show()

if __name__ == '__main__':

    anim = True
    dt = 0.025
    file = 'abc.dat'
    interp = False
    xmax = 2.5

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ad:f:ix:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            anim = not anim
        elif o == "-d":
            dt = float(a)
        elif o == "-f":
            file = a
        elif o == "-i":
            interp = not interp
        elif o == "-x":
            xmax = float(a)
        else:
            print "unexpected argument", o
            sys.exit(1)

    nt, t, i, x, y, z = read_data(file)

    # Interpolate to the desired dt.

    if interp:
        nt, t, x, y = interpolate_data(t, x, y, dt)

    # Animate or plot the data.

    if anim:
        animate_data(t, x, y, i[:nt], 'x', 'y', xmax)
    else:
        plot_data(t, x, 't', 'x')
