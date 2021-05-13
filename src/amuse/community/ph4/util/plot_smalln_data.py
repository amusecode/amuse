import sys
import getopt
import math
import numpy
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animate

# TODO:	    Colormap should respond to changes in index array.

def unpack_line(r, np, nq):
    
    tt = float(r[0])
    ii = []
    mm = []
    xx = []
    yy = []
    zz = []
    for j in range(np):
        ii.append(int(r[nq*j+1]))
        if nq == 4:
            mm.append(2.**(-j))
        else:
            mm.append(float(r[nq*j+2]))
        xx.append(float(r[nq*j+nq-2]))
        yy.append(float(r[nq*j+nq-1]))
        zz.append(float(r[nq*j+nq]))

    return tt, ii, mm, xx, yy, zz

def interpolate(xp, xx, tfac):
    
    xint = []
    for j in range(len(xx)):
        xint.append(xp[j] + tfac*(xx[j]-xp[j]))
    return xint

def read_data(filename, dt, N, have_mass):

    # Read in the data file.  Format is defined in smalln.

    print('reading data from', filename)
    sys.stdout.flush()

    f = open(filename)
    ll = f.readlines()
    f.close()

    nt = len(ll)
    print(nt, 'records read')
    sys.stdout.flush()

    # Do some preliminary analysis.

    r = ll[0].split()			# trust the first line
    if have_mass:
        nq = 5
    else:
        nq = 4
    np = (len(r)-1)/nq
    print('np =', np, 'nq =', nq)
    
    tt, ii, mm, xx, yy, zz = unpack_line(r, np, nq)
    tnext = tt

    if N > 0:				# N trumps dt
        r = ll[nt-1].split()		# last line
        tlast, dum, dum, dum, dum = unpack_line(r, np, nq)
        pdt = False
        if dt > 0: pdt = True
        dt = (tlast-tt)/N
        if pdt: print('resetting dt =', dt)

    t = []
    i = []
    m = []
    x = []
    y = []
    z = []
    line = 0
    short = 0

    for l in ll:

        r = l.split()
        line += 1

	# Format is: time  id1 m1 x1 y1 z1  id2 m2 x2 y2 z2  id3 m3 x3 y3 z3 ...

        if nq*np + 1 == len(r):

            # Save old data and unpack the line.

            tp = tt
            mp= mm
            xp = list(xx)
            yp = list(yy)
            zp = list(zz)
            tt, ii, mm, xx, yy, zz = unpack_line(r, np, nq)

            while tt >= tnext:

                if line == 1 or dt == 0.0:

                    tint = tt
                    mint = mm
                    xint = list(xx)
                    yint = list(yy)
                    zint = list(zz)
                    
                else:
                    
                    # Interpolate to time tnext.

                    tint = tnext
                    tfac = (tnext-tp)/(tt-tp)
                    mint = interpolate(mp, mm, tfac)
                    xint = interpolate(xp, xx, tfac)
                    yint = interpolate(yp, yy, tfac)
                    zint = interpolate(zp, zz, tfac)

                # Append the data to the master arrays.

                t.append(tint)
                i.append(ii)
                m.append(mint)
                x.append(xint)
                y.append(yint)
                z.append(zint)

                tnext += dt
                if dt == 0.0: break

        else:				# temporary 2-body treatment
            
            #print 'ignoring short record:', len(r), '!=', 4*np+1
            short += 1

    ta = numpy.array(t)
    ia = numpy.array(i)
    ma = numpy.array(m)
    xa = numpy.array(x)			# should be nt x np
    ya = numpy.array(y)
    za = numpy.array(z)

    nt = len(ta)
    print('nt =', nt)
    print(short, 'short records')
    #print 'xa.shape =', xa.shape	# should be nt x np
    return nt, ta, ia, ma, xa, ya, za

def print_help():
    print('keyboard controls:')
    print('    space    pause/resume')
    print('    a        expand view to contain all particl;es')
    print('    h        print this message')
    print('    q        quit')
    print('    z        zoom in')
    print('    Z        zoom out')
    print('    right    pan right')
    print('    left     pan left')
    print('    up       pan up')
    print('    down     pan down')
    print('    <        first frame')
    print('    >        last frame')
    print('mouse click reverses direction')

# We seem to need some of these global for the animate functions to work...

nf = 0
df = 1
xmin = 0.0
xmax = 0.0
ymin = 0.0
ymax = 0.0
anim_running = True
current_frame = 0
shift = 0.25
zoom = 1.5

def animate_data(t, m, x, y, id, lx, ly, scale, delay):

    global xmin, xmax, ymin, ymax

    print('animating data')
    colormap = ['r', 'y', 'b', 'm', 'g', 'c', 'k', 'k', 'k', 'k', 'k', 'k']

    # Determine length scales.

    xmin = numpy.min(x[0,:])
    xmax = numpy.max(x[0,:])
    dx = xmax - xmin
    xav = 0.5*(xmin+xmax)
    ymin = numpy.min(y[0,:])
    ymax = numpy.max(y[0,:])
    dy = ymax - ymin
    yav = 0.5*(ymin+ymax)

    # Allow modification of overall scale by user-specified scale factor.

    dx = 5*max(dx, dy)*scale
    xmin = xav - 0.5*dx
    xmax = xav + 0.5*dx
    ymin = yav - 0.5*dx
    ymax = yav + 0.5*dx

    # Determine time scales.

    dtscale = t[-1]
    if t[0] < 0: dtscale = max(dtscale, -t[0])
    nf = int(math.floor(math.log10(dtscale)))
    #print 'dtscale =', dtscale, 'nf =', nf
    if nf > 3:
        tformat = '%'+str(nf+2)+'.0f'
    elif nf > 0:
        tformat = '%'+str(nf+4)+'.2f'
    else:
        tformat = '%'+str(-nf+7)+'.'+str(-nf+4)+'f'
    lformat = tformat+'/'+tformat+'    frame %d/%d'
    s = numpy.log(list(m[0,:]))
    smin = numpy.min(s)
    smax = numpy.max(s)
    s = 5 + 25*(s-smin)/(smax-smin)    # sizes logarithmic in mass, range 5-30
    print('s =', s)
    nt,np = x.shape
    fig = plt.figure()
    
    scat = plt.scatter(x[0,:], y[0,:], c=colormap[:np], s=s)
    #plt.axis('square')	# fails in ubuntu
    plt.xlabel(lx)
    plt.ylabel(ly)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.grid()

    def update(frame):

        # What idiot decided that the syntax for set_offsets should
        # be so different from the syntax for scatter??

        off = []
        for j in range(np):
            off.append(x[frame,j])
            off.append(y[frame,j])
        scat.set_offsets(off)
        text = fig.suptitle('')
        text.set_text('')
        text.set_text('time '+lformat%(t[frame],t[-1], frame, nt-1))

        return scat,

    def nextframe(nt):			# step nf forward or backward by df
        global nf, df			# stay fixed at nf = 0 or nt-1
        global current_frame
        nf = 0
        while nf < nt:
            current_frame = nf
            yield nf
            nf += df
            if nf < 0: nf = 0
            if nf >= nt: nf = nt-1
    
    def new_limits(x, y, fac):
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        xm = 0.5*(xmin+xmax)
        dx = xmax - xm
        xmax = xm + fac*dx
        xmin = xm - fac*dx
        ym = 0.5*(ymin+ymax)
        dy = ymax - ym
        ymax = ym + fac*dy
        ymin = ym - fac*dy
        return xmin, xmax, ymin, ymax

    def onClick(event):			# reverse direction on mouse click
        global nf, df
        df = -df

    def onKey(event):			# manage key presses
        global anim_running
        global nf, df
        global xmin, xmax, ymin, ymax
        global current_frame
        global shift, zoom
        if event.key == 'a':		# a = reset limits to current particles
            xmin, xmax, ymin, ymax = \
			new_limits(x[current_frame],
                                   y[current_frame], 1.25)
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
        elif event.key == 'q':		# q = quit
            sys.exit(0)
        elif event.key == 'z':		# z = zoom in
            xm = 0.5*(xmin+xmax)
            dx = xmax - xm
            xmax = xm + dx/zoom
            xmin = xm - dx/zoom
            ym = 0.5*(ymin+ymax)
            dy = ymax - ym
            ymax = ym + dy/zoom
            ymin = ym - dy/zoom
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
        elif event.key == 'Z':		# Z = zoom out
            xm = 0.5*(xmin+xmax)
            dx = xmax - xm
            xmax = xm + dx*zoom
            xmin = xm - dx*zoom
            ym = 0.5*(ymin+ymax)
            dy = ymax - ym
            ymax = ym + dy*zoom
            ymin = ym - dy*zoom
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
        elif event.key == ' ':		# space = pause/restart
            if anim_running:
                anim.event_source.stop()
                anim_running = False
            else:
                anim.event_source.start()
                anim_running = True
        elif event.key == 'right':
            dx = shift*(xmax - xmin)
            xmax += dx
            xmin += dx
            plt.xlim(xmin, xmax)
        elif event.key == 'left':
            dx = shift*(xmax - xmin)
            xmax -= dx
            xmin -= dx
            plt.xlim(xmin, xmax)
        elif event.key == 'up':
            shift = 0.5
            dy = shift*(ymax - ymin)
            ymax += dy
            ymin += dy
            plt.ylim(ymin, ymax)
        elif event.key == 'down':
            shift = 0.5
            dy = shift*(ymax - ymin)
            ymax -= dy
            ymin -= dy
            plt.ylim(ymin, ymax)
        elif event.key == '<':
            nf = 0
            df = 1
        elif event.key == '>':
            nf = nt - 1
            df = -1
        elif event.key == 'h':
            print_help()
        else:
            print('key =', event.key)

    fig.canvas.mpl_connect('key_press_event', onKey)
    fig.canvas.mpl_connect('button_press_event', onClick)
    anim = animate.FuncAnimation(fig, update, frames=nextframe(nt),
                                 interval=delay, repeat=False)
    plt.show()

def plot_data(t, x, lx, ly):
    nt,np = x.shape
    fig = plt.figure()
    plt.plot(t, x, 'b-')
    plt.xlabel(lx)
    plt.ylabel(ly)
    plt.grid()
    plt.show()

def write_data(t, x, y, i, filename):
    nt,np = x.shape
    f = open(filename)
    for i in range(nt):
        s = str(t[i])+' '
        for j in range(np):
            s += str(i(i,j))+' '+str(x(i,j))+' '+str(y(i,j))+' 0.0\n'
        f.write(s)
    f.close()

if __name__ == '__main__':

    anim = True
    dt = 0.0
    delay = 30
    file = 'abc.dat'
    have_mass = True
    N = 0
    outfile = None
    proj = 3
    scale = 1.5

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ad:D:f:mN:o:p:s:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            anim = not anim
        elif o == "-d":
            dt = float(a)
        elif o == "-D":
            delay = int(a)
        elif o == "-f":
            file = a
        elif o == "-m":
            have_mass = not have_mass
        elif o == '-N':
            N = int(a)
        elif o == "-o":
            outfile = a
        elif o == '-p':
            proj = int(a)
        elif o == "-s":
            scale = float(a)
        else:
            print("unexpected argument", o)
            sys.exit(1)

    nt, t, i, m, x, y, z = read_data(file, dt, N, have_mass)

    # Animate or plot the data.

    if anim:
        a1 = x
        a2 = y
        l1 = 'x'
        l2 = 'y'
        if proj == 1:
            a1 = y
            a2 = z
            l1 = 'y'
            l2 = 'z'
        elif proj == 2:
            a2 = z
            l2 = 'z'
        animate_data(t, m, a1, a2, i[:nt], l1, l2, scale, delay)
    else:
        plot_data(t, m, x, 't', 'x')
